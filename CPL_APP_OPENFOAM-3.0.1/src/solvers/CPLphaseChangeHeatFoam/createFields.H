Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField delU
(
    IOobject
    (
        "delU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;

    autoPtr<phaseChangeTwoPhaseMixture> twoPhaseProperties =
        phaseChangeTwoPhaseMixture::New(U, phi);

    	const dimensionedScalar& rho1 = twoPhaseProperties->rho1();
        const dimensionedScalar& rho2 = twoPhaseProperties->rho2();

        const dimensionedScalar& K1 = twoPhaseProperties->K1();
        const dimensionedScalar& K2 = twoPhaseProperties->K2();

        const dimensionedScalar& C1 = twoPhaseProperties->C1();
        const dimensionedScalar& C2 = twoPhaseProperties->C2();

        const volScalarField& TSatLocal = twoPhaseProperties->TSatLocal();
        const dimensionedScalar& pSat = twoPhaseProperties->pSat();



immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
/*
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
*/
const dimensionedScalar& neu1 = mixture.neu1();
const dimensionedScalar& neu2 = mixture.neu2();


// Need to store rho for ddt(rho, U)

/*
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField alpha2
    (
        IOobject
        (
            "alpha2",
            runTime.timeName(),
            mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh
    );
*/
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + (1-alpha1)*rho2,
    alpha1.boundaryField().types()
);
rho.oldTime();

volScalarField neu
(
    IOobject
    (
        "neu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    (alpha1*neu1*rho1 + (1-alpha1)*neu2*rho2)/(alpha1*rho1 + (1-alpha1)*rho2),
    alpha1.boundaryField().types()
);

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);



// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);




#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

/*volScalarField rgh
(
    IOobject
    (
        "rgh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*gh
);*/

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

#include "readTransportProperties.H"

    Info<< "Creating field kinetic energy K\n" << endl;
        volScalarField K("K", 0.5*magSqr(U));

        volScalarField mDotT
           (
           		IOobject
                   (
                       "mDotT",
                       runTime.timeName(),
                       mesh,
                       IOobject::NO_READ,
                       IOobject::AUTO_WRITE
                   ),
                mesh,
                dimensionedScalar ("mDotT",dimDensity/dimTime,0)
           );

        volScalarField mDotP
                  (
                  		IOobject
                          (
                              "mDotP",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ,
                              IOobject::AUTO_WRITE
                          ),
                       mesh,
                       dimensionedScalar ("mDotP",dimDensity/dimTime,0)
                  );
        volScalarField mDotAlpha
                  (
                  		IOobject
                          (
                              "mDotAlpha",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ,
                              IOobject::AUTO_WRITE
                          ),
                       mesh,
                       dimensionedScalar ("mDotAlpha",dimDensity/dimTime,0)
                  );
// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;


