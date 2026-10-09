/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    CPLInterFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Modified from interFoam to couple it with LAMMPS solver using CPL Library

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "CorrectPhi.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "CPLSocketFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Check if coupled based on cpl/COUPLER.in input file
    bool coupled;
    if (file_exists("./cpl/COUPLER.in")) {
        Info<< "Assuming coupled run as cpl/COUPLER.in exists\n" << endl;
        coupled=true;
    } else {
        Info<< "Assuming uncoupled run as cpl/COUPLER.in does not exist\n" << endl;
        coupled=false;
    }


    // Create a CPL object (not used if uncoupled) and intialises MPI
    CPLSocketFOAM CPL;
    if (coupled)
        CPL.initComms(argc, argv);

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // MPI_Init is called somewhere in the PStream library
    if (coupled)
        CPL.init2pCFD(runTime, mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	

	// Initial communication to initialize domains
    if (coupled){
        CPL.pack(U, p, alpha1, neu, mesh, CPL.VEL | CPL.ALPHA);
        CPL.send();
        CPL.recvTwoPhase();
        CPL.unpackTwoPhase(U, delU, alpha1, rho1, rho2, mesh);
    //    CPL.recvVelocityPressure();
    //    CPL.unpackVelocityPressure(U, p, mesh);
    }

    while (runTime.run())
    {
    
        if (coupled){
        CPL.pack(U, p, alpha1, neu, mesh, CPL.VEL | CPL.ALPHA);
        CPL.send();
        CPL.recvTwoPhase();
        CPL.unpackTwoPhase(U, delU, alpha1, rho1, rho2, mesh);
    //    CPL.recvVelocityPressure();
    //    CPL.unpackVelocityPressure(U, p, mesh);
    }
    
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }
        


        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    CPL::finalize();

    return 0;
}


// ************************************************************************* //
