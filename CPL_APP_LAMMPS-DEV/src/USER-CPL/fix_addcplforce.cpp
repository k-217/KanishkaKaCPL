/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <string>
#include "fix_addcplforce.h"
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "variable.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddCPLForce::FixAddCPLForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), estr(NULL), idregion(NULL), sforce(NULL)

{
  if (narg < 7) error->all(FLERR,"Illegal fix addcplforce command 1");
  int nlocal = atom->nlocal;
  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  virial_flag = 1;

  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }
  int flag;
  int indexH = atom->find_custom("hRef", flag);
    	double *dvectorH = atom->dvector[indexH];
  for (int i = 0; i < nlocal; ++i)
    {
    	dvectorH[i] = force->numeric(FLERR,arg[6]);
    	}

  bndryAvgOld = 0;
//std::cout << "href " << std::endl;
//    href = new double[2];
//    href[0] = force->numeric(FLERR,arg[6]);
  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addcplforce command 2");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix addcplforce command 3");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addcplforce command 4");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix addcplforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addcplforce command 5");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix addcplforce command 6");
      iarg += 2;
    } 
    //std::cout << "end of constructor 1" << std::endl;
    else error->all(FLERR,"Illegal fix addcplforce command 7");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"addcplforce:sforce");
  std::cout << "end of constructor" << std::endl;
}

/* ---------------------------------------------------------------------- */

FixAddCPLForce::~FixAddCPLForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
//  delete [] href;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixAddCPLForce::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddCPLForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix addcplforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix addcplforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix addcplforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix addcplforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix addcplforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix addcplforce is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix addcplforce does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix addcplforce is invalid style");
  } else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix addcplforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant force in fix addcplforce");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix addcplforce");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

   std::ofstream foutput;
   foutput.open ("href.txt",std::ios::app); 
 
   foutput<< "delHref " << "bndryAvg " <<"href" << std::endl;
   foutput.close(); 
}

/* ---------------------------------------------------------------------- */

void FixAddCPLForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddCPLForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddCPLForce::post_force(int vflag)
{
  double **x = atom->x;
  double **vel = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6],delHref;
  double bndryAvg=0;
  int nlocal = atom->nlocal;
  int flag;
  int ncount=0;
  int index = atom->find_custom("alpha", flag);
  double *dvector = atom->dvector[index];


	
  int indexH = atom->find_custom("hRef", flag);
  double *dvectorH = atom->dvector[indexH];

  
  if (update->ntimestep % nevery) return;
	if (!(update->ntimestep % 500)) std::cout << "href = " << dvectorH[0] << std::endl;
  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if ((varflag == ATOM || estyle == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"addcplforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
        foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2];
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        double bndryForce=0,bndryX, bndryZ, x0, x1, x2, cutH;
        
        
        
//        if(update->ntimestep<1000)
//	dvector[i] = 1;
	
	
//        std::cout << "href = " << href << std::endl;
/*        if(x[i][1]-dvectorH[i] <3.5)
          {
          bndryForce = (0.4078285691*pow((x[i][1]-dvectorH[i]),2)+4.2849123468*(x[i][1]-dvectorH[i])+11.3858621732);
          }
          else
          {
          bndryForce = (-0.0065655868*pow((x[i][1]-dvectorH[i]),3)+0.0366196158*pow((x[i][1]-dvectorH[i]),2)-0.1154670851*(x[i][1]-dvectorH[i])+0.1736053968);
          }*/
	/*if(dvectorH[i]-x[i][1] <0)
	  {
		bndryForce = 0;
		}
        else if(dvectorH[i]-x[i][1] <1.7)
          {
          bndryForce = pow(10,-10)*(0.3827489048*(dvectorH[i]-x[i][1])+0.1976351964);
          }
          else if(dvectorH[i]-x[i][1] <8)
          {
          bndryForce = pow(10,-10)*(-0.001212422*pow((dvectorH[i]-x[i][1]),3)+0.0397551787*pow((dvectorH[i]-x[i][1]),2)-0.4270250996*(dvectorH[i]-x[i][1])+1.497557021);
          } 
          else if(dvectorH[i]-x[i][1] <8)
          {
		bndryForce = 0;
		}*/


/*	if(dvectorH[i]-x[i][1] <-0.625)
	  {
		bndryForce = 0;
		}
        else if(dvectorH[i]-x[i][1] <-0.125)
          {
          bndryForce =0;//(1.8400738*(dvectorH[i]-x[i][1])+1.2187945417);//(7.4766388*(dvectorH[i]-x[i][1])+4.8994164167);
          }
          else if(dvectorH[i]-x[i][1] <0.5)
          {
          bndryForce =-0.5;//(-2.899924*(dvectorH[i]-x[i][1])+0.5781519583);// (-8.53649*(dvectorH[i]-x[i][1])+2.4402415833);
          }
          else if(dvectorH[i]-x[i][1] <2.875)
          {
          bndryForce =(0.2956725781*pow((dvectorH[i]-x[i][1]),3)-1.4994818713*pow((dvectorH[i]-x[i][1]),2)+2.7254065281*(dvectorH[i]-x[i][1])-1.9228269052);
          }
          else if(dvectorH[i]-x[i][1] <7)
          {
          bndryForce =(-0.0590644555*pow((dvectorH[i]-x[i][1]),3)+0.9521026206*pow((dvectorH[i]-x[i][1]),2)-4.9396617598*(dvectorH[i]-x[i][1])+8.1884596047);
          }
          else
          {
		bndryForce = 0;
		}*/








		
	/*if(dvectorH[i]-x[i][1] <0)
	  {
		bndryForce = -5;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
	if((dvectorH[i]-x[i][1]) <3)
	  {
		bndryForce = (-1.1055763868*(dvectorH[i]-x[i][1])+2.4698660286);
		
		}
	else if(abs(dvectorH[i]-x[i][1]) <6)
          {
          bndryForce = (-0.06409966613*pow((dvectorH[i]-x[i][1]),2)+0.8323267514*(dvectorH[i]-x[i][1])-2.7158572434);
          }
        else
          {
          bndryForce =0;
          } 
*/



//OPPOSITE FORCE
//std::cout << "href = " << href << std::endl;
/*	if(dvector[i]>0.8 && abs(dvectorH[i]-x[i][1]) <6)
	{
		if(abs(dvectorH[i]-x[i][1]) <3)
		  {
			bndryForce = -(-1.1055763868*abs(dvectorH[i]-x[i][1])+2.4698660286);
		
			}
		else
		  {
			bndryForce = -(-0.06409966613*pow((dvectorH[i]-x[i][1]),2)+0.8323267514*abs(dvectorH[i]-x[i][1])-2.7158572434);
			}
          }
        else
          {
          bndryForce =0;
          } */
//	std::cout << "alpha = " << dvector[i] << std::endl;





//alpha based model




/*	cutH = 4.51875*dvector[i];
	if((dvector[i]<0.1))
	bndryForce = 0;
	else
	{
	if((dvectorH[i]-x[i][1])<2.875-cutH)
	  {
		x2=0.3041886672*pow(cutH,2)-0.1223429689*cutH+0.0659099843;
		x1=0.5742383391*pow(cutH,2)-1.5903295649*cutH+0.2656891561;
		x0=0.066092740*pow(cutH,3)+0.2072779714*pow(cutH,2)-0.120223294*cutH+0.0292459177;
		bndryForce = -(x2*pow((dvectorH[i]-x[i][1]),2)+x1*(dvectorH[i]-x[i][1])+x0);
		}
	else if((dvectorH[i]-x[i][1])>2.3925)
	  {
		if(dvector[i]<0.7)
		{
			x2=0.0237964269*pow(cutH,3)-0.1939125369*pow(cutH,2)+0.5306296611*cutH+0.065911743;
			x1=-0.1442490485*pow(cutH,3)+1.1678472144*pow(cutH,2)-3.1740085729*cutH-0.494253892;
			x0=0.20369724*pow(cutH,3)-1.58358650221*pow(cutH,2)+4.0631975186*cutH+0.878992488;
			bndryForce = -(x2*pow((dvectorH[i]-x[i][1]),2)+x1*(dvectorH[i]-x[i][1])+x0);
			}
		else
		bndryForce = -(0.5871238095*pow((dvectorH[i]-x[i][1]),2)-3.5961529524*(dvectorH[i]-x[i][1])+4.6075690714);
		if(bndryForce < 0)
		bndryForce = 0;
		}
	else
	  {
		if(dvector[i]<0.7)
		{
			x2=0.2314569562*pow(cutH,2)-0.9293813532*cutH-0.1562357024;
			x1=-1.0291059573*pow(cutH,2)+2.5898688502*cutH+1.9244996471;
			x0=0.8990125367*pow(cutH,2)-0.7798248997*cutH-3.9117006279;
			bndryForce = -(x2*pow((dvectorH[i]-x[i][1]),2)+x1*(dvectorH[i]-x[i][1])+x0);
			}
		else
		bndryForce = -(-1.3065624867*(dvectorH[i]-x[i][1])+2.3848183419);
		}





	}*/




	if(dvector[i]>0.8)
	{
		if((dvectorH[i]-x[i][1])>2.3925)
			{
			bndryForce = -(0.5871238095*pow((dvectorH[i]-x[i][1]),2)-3.5961529524*(dvectorH[i]-x[i][1])+4.6075690714);
			if(bndryForce < 0)
				bndryForce = 0;
		}
		else
		bndryForce = -(-1.3065624867*(dvectorH[i]-x[i][1])+2.3848183419);
		}
	else
	bndryForce = 0;

//	bndryForce = bndryForce*dvector[i];


//High density 0.03665
/*	if(dvector[i]>0.55)
		{
		if(abs(dvectorH[i]-x[i][1]) <2)
		  {
			if((dvector[i]<0.75 && dvectorH[i]-x[i][1]) <0.75)
			bndryForce = -1.5;
			else
			bndryForce = -(-1.46052461904762*abs(dvectorH[i]-x[i][1])+1.58958074404762);
			}
		else if(abs(dvectorH[i]-x[i][1]) <4)
		  {
			  bndryForce = -(-0.58989245021645*pow((dvectorH[i]-x[i][1]),2)+3.15611805541125*abs(dvectorH[i]-x[i][1])-2.74879766011904);
			  }
		else if(abs(dvectorH[i]-x[i][1]) <6)
		  {
			  bndryForce = -(-0.064126496969696*pow((dvectorH[i]-x[i][1]),2)+0.930099190303023*abs(dvectorH[i]-x[i][1])-3.31979818227271);
			  }
		else
		  {
			  bndryForce =0;
			  }
		}
	else if(dvector[i]>0.3)
	  {
		if(abs(dvectorH[i]-x[i][1]) <0.625)
		  {
			bndryForce = 0.2;
			}
		else if(abs(dvectorH[i]-x[i][1]) <2.125)
		  {
			  bndryForce = -(-1.0245918095*pow((dvectorH[i]-x[i][1]),2)+2.701990476*abs(dvectorH[i]-x[i][1])-1.4557636875);
			  }
		else if(abs(dvectorH[i]-x[i][1]) <5)
		  {
			  bndryForce = -(0.2606535758*pow((dvectorH[i]-x[i][1]),2)-1.6618800727*abs(dvectorH[i]-x[i][1])-1.931717152);
			  }
		else
		  {
			  bndryForce =0;
			  }
		}
	else if(dvector[i]>0.1)
	  {
		if(abs(dvectorH[i]-x[i][1]) <2)
		  {
			  bndryForce = -(0.3180209004*pow((dvectorH[i]-x[i][1]),2)-0.8220317593*abs(dvectorH[i]-x[i][1])+0.020907263);
			  }
		else if(abs(dvectorH[i]-x[i][1]) <4)
		  {
			  bndryForce = -(0.3450595238*pow((dvectorH[i]-x[i][1]),2)-2.1703345714*abs(dvectorH[i]-x[i][1])+2.8857561057);
			  }
		else
		  {
			  bndryForce =0;
			  }
		}
	else 
	  {
	  bndryForce =0;
			  
		}*/






/*	if(dvectorH[i]-x[i][1] <1)
	  {
		bndryForce =0.05*(2*(dvectorH[i]-x[i][1])); //0.05*(-1.1055763868*(dvectorH[i]-x[i][1])+2.4698660286);
		
		}
	else if(dvectorH[i]-x[i][1] <2.5)
          {
          bndryForce = 0.05*(-2*(dvectorH[i]-x[i][1])/3+4); //0.05*(-0.06409966613*pow((dvectorH[i]-x[i][1]),2)+0.8323267514*(dvectorH[i]-x[i][1])-2.7158572434);
          }
        else if(dvectorH[i]-x[i][1] <6)
          {
          bndryForce =0.05*((dvectorH[i]-x[i][1])/3.5-6/3.5);
          }
	else
          {
          bndryForce =0;
          }
*/



//          std::cout << "bndryForce = " << bndryForce << std::endl;
          //if(dvector[i]==1)	std::cout << "Recv dvector[" << i << "] = " << dvector[i] << " x = " << x[i][0] << " z = " << x[i][2] << std::endl;
	/*if(dvector[i]<0.99) 
	{
		if(dvectorH[i]-x[i][1] <0)
	 	{
		bndryForce = -2;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
		if((dvectorH[i]-x[i][1]) <2.625 && (dvectorH[i]-x[i][1])>0)
		 bndryForce= ((0.0244*pow((dvectorH[i]-x[i][1]),2))+0.0013*(dvectorH[i]-x[i][1])-0.1946);
		else
		 bndryForce =0;
	}
	else
	{
		if(dvectorH[i]-x[i][1] <0)
	  {
		bndryForce = -5;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
	if((dvectorH[i]-x[i][1]) <3)
	  {
		bndryForce = (-1.1055763868*(dvectorH[i]-x[i][1])+2.4698660286);
		
		}
	else if(abs(dvectorH[i]-x[i][1]) <6)
          {
          bndryForce = (-0.06409966613*pow((dvectorH[i]-x[i][1]),2)+0.8323267514*(dvectorH[i]-x[i][1])-2.7158572434);
          }
        else
          {
          bndryForce =0;
          } 
	} */



//          std::cout << "bndryForce = " << bndryForce << std::endl;
          //if(dvector[i]==1)	std::cout << "Recv dvector[" << i << "] = " << dvector[i] << " x = " << x[i][0] << " z = " << x[i][2] << std::endl;






//Argon equation

/*	if(dvector[i]<0.99) 
	{
//		if(dvectorH[i]-x[i][1] <0)
	 	{
		bndryForce = -0.23;
//		vel[i][1] = -vel[i][1];
//		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
//		else 
		if(abs(dvectorH[i]-x[i][1]) <3.5)
		 bndryForce= -(-0.0401343547*pow((dvectorH[i]-x[i][1]),2)-0.0028221134*(dvectorH[i]-x[i][1])+0.2281826201);
		else if(abs(dvectorH[i]-x[i][1]) <6.2)
		 bndryForce =-(-0.0260461073*pow((dvectorH[i]-x[i][1]),2)+0.3302130055*(dvectorH[i]-x[i][1])-1.0987189292);
		else if(abs(dvectorH[i]-x[i][1]) <9.0)
		 bndryForce =-(-0.0044658909*pow((dvectorH[i]-x[i][1]),2)+0.0797119557*(dvectorH[i]-x[i][1])-0.3674322064);
		else if(abs(dvectorH[i]-x[i][1]) <15.0)
		 bndryForce= 0;//-(0.001336637*(dvectorH[i]-x[i][1])-0.0200531781);
		else 
		 bndryForce= 0;
	}
	else
	{
		bndryForce = 0;
		if(dvectorH[i]-x[i][1] <0)
	  {
		bndryForce = -5;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
	} 

	bndryForce = 0.0221930911*bndryForce/0.0287;

*/


//Nitrogen equation

/*	if(dvector[i]<0.99) 
	{
		if(dvectorH[i]-x[i][1] <0)
	 	{
		bndryForce = -6;
//		vel[i][1] = -vel[i][1];
//		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
		else if((dvectorH[i]-x[i][1]) <2.8)
		 bndryForce= -0.01*((0.7175*pow((dvectorH[i]-x[i][1]),2))-0.0844*(dvectorH[i]-x[i][1])+5.8956);
		else if((dvectorH[i]-x[i][1]) <3.0)
		 bndryForce =-(-0.9588*(dvectorH[i]-x[i][1])+2.808);
		else if((dvectorH[i]-x[i][1]) <6.0)
		 bndryForce= -((-0.0068*pow((dvectorH[i]-x[i][1]),2))+0.0883*(dvectorH[i]-x[i][1])-0.2968);
		else if((dvectorH[i]-x[i][1]) <10.0)
		 bndryForce= -((-0.0005*pow((dvectorH[i]-x[i][1]),2))+0.0103*(dvectorH[i]-x[i][1])-0.0538);
		else 
		 bndryForce= 0;
	}
	else
	{
		if(dvectorH[i]-x[i][1] <0)
	  {
		bndryForce = -5;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
	if((dvectorH[i]-x[i][1]) <3)
	  {
		bndryForce = (-1.1055763868*(dvectorH[i]-x[i][1])+2.4698660286);
		
		}
	else if(abs(dvectorH[i]-x[i][1]) <6)
          {
          bndryForce = (-0.06409966613*pow((dvectorH[i]-x[i][1]),2)+0.8323267514*(dvectorH[i]-x[i][1])-2.7158572434);
          }
        else
          {
          bndryForce =0;
          } 
	} */


//Two-phase bubble eqution

/*
	if(dvector[i]<0.99) 
	{
		if(dvectorH[i]-x[i][1] <0)
	 	{
		bndryForce = -2;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
		if((dvectorH[i]-x[i][1]) <2.625 && (dvectorH[i]-x[i][1])>0)
		 bndryForce=-1 * ((0.0244*pow((dvectorH[i]-x[i][1]),2))+0.0013*(dvectorH[i]-x[i][1])-0.1946);
		else
		 bndryForce =0;
	}
	else
	{
		if(dvectorH[i]-x[i][1] <0)
	  {
		bndryForce = -5;
		vel[i][1] = -vel[i][1];
		x[i][1] = 2*dvectorH[i]-x[i][1];
		}
	if((dvectorH[i]-x[i][1]) <3)
	  {
		bndryForce = (-1.1055763868*(dvectorH[i]-x[i][1])+2.4698660286);
		
		}
	else if(abs(dvectorH[i]-x[i][1]) <6)
          {
          bndryForce = (-0.06409966613*pow((dvectorH[i]-x[i][1]),2)+0.8323267514*(dvectorH[i]-x[i][1])-2.7158572434);
          }
        else
          {
          bndryForce =0;
          } 
	}
	
*/

        f[i][0] += (xvalue);
        f[i][1] += (yvalue+bndryForce);
        f[i][2] += (zvalue);


	if(bndryForce!=0 && dvector[i] > 0.9)
	{
	
	ncount++;
	bndryAvg += (bndryForce-bndryAvg)/ncount;
	}

        if (evflag) {
          v[0] = xvalue * unwrap[0];
          v[1] = yvalue * unwrap[1];
          v[2] = zvalue * unwrap[2];
          v[3] = xvalue * unwrap[1];
          v[4] = xvalue * unwrap[2];
          v[5] = yvalue * unwrap[2];
          v_tally(i,v);
        }
        
      }

  // variable force, wrap with clear/add
  // potential energy = evar if defined, else 0.0
  // wrap with clear/add

  } else {
    double unwrap[3];

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);
    if (estyle == ATOM)
      input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
        if (xstyle == ATOM) xvalue = sforce[i][0];
        if (ystyle == ATOM) yvalue = sforce[i][1];
        if (zstyle == ATOM) zvalue = sforce[i][2];

        if (estyle == ATOM) {
          foriginal[0] += sforce[i][3];
        } else {
          if (xstyle) foriginal[0] -= xvalue*unwrap[0];
          if (ystyle) foriginal[0] -= yvalue*unwrap[1];
          if (zstyle) foriginal[0] -= zvalue*unwrap[2];
        }
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];

        if (xstyle) f[i][0] += xvalue;
        if (ystyle) f[i][1] += yvalue;
        if (zstyle) f[i][2] += zvalue;
        if (evflag) {
          v[0] = xstyle ? xvalue*unwrap[0] : 0.0;
          v[1] = ystyle ? yvalue*unwrap[1] : 0.0;
          v[2] = zstyle ? zvalue*unwrap[2] : 0.0;
          v[3] = xstyle ? xvalue*unwrap[1] : 0.0;
          v[4] = xstyle ? xvalue*unwrap[2] : 0.0;
          v[5] = ystyle ? yvalue*unwrap[2] : 0.0;
          v_tally(i, v);
        }
      }
    }
  }

	bndryAvg=bndryAvg*0.1+0.9*bndryAvgOld;
	bndryAvgOld = bndryAvg;
	delHref = -3.2251513876*bndryAvg;//+0.0941031303;//-0.9738093081*bndryAvg*bndryAvg-2.78788531*bndryAvg-0.1170954947;//-0.05
	

//-1.9157317931*bndryAvg+0.3809782621;//-14.9710304146*bndryAvg-4.964500164;//-25.5190923422*bndryAvg-1.7483391969;

    for (int i = 0; i < nlocal; i++) {
//      if (mask[i] & groupbit) 
	{dvectorH[i]+=delHref;
      }
    }

 std::ofstream foutput;
 foutput.open ("href.txt",std::ios::app); 
 
   foutput<< delHref << " " << bndryAvg << " " << dvectorH[0] << " " << ncount << std::endl;
   foutput.close(); 

//   std::cout << " delHref " <<  delHref << " bndryAvg " << bndryAvg << " dvectorH[i]" << dvectorH[0] << std::endl;
}

/* ---------------------------------------------------------------------- */

void FixAddCPLForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddCPLForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddCPLForce::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddCPLForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAddCPLForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
