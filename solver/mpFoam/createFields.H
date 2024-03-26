#include "createRDeltaT.H"

Info<< "Reading personlized input parameters" << endl;
//IOdictionary solid1Inputs
//(
//    IOobject
//    (
//        "solid1Inputs",
//        runTime.constant(),
//        mesh,
//        IOobject::MUST_READ_IF_MODIFIED,
//        IOobject::NO_WRITE
//    )
//);

IOdictionary solid1Inputs
(
    IOobject
    (
        "solid1Inputs",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary solid2Inputs
(
    IOobject
    (
        "solid2Inputs",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar molW1("Mv", dimMass/dimMoles, solid1Inputs);
const dimensionedScalar Csat1("Csat", dimMoles/dimVolume, solid1Inputs);
const dimensionedScalar Cinit1("Cinit", dimMoles/dimVolume, solid1Inputs);
const word precDissSolid1(solid1Inputs.lookup("precOrDiss"));
const dimensionedScalar molW2("Mv", dimMass/dimMoles, solid2Inputs);
const dimensionedScalar Csat2("Csat", dimMoles/dimVolume, solid2Inputs);
const dimensionedScalar Cinit2("Cinit", dimMoles/dimVolume, solid2Inputs);
const word precDissSolid2(solid2Inputs.lookup("precOrDiss"));
const bool nucleationFlag1(solid1Inputs.lookupOrDefault<bool>("nucleation", false));

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

#include "createPhi.H"

Info<< "Reading field C1\n" << endl;
volScalarField C1
(
    IOobject
    (
        "C1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C2\n" << endl;
volScalarField C2
(
    IOobject
    (
        "C2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating reaction source terms for C and alpha\n" << endl;
volScalarField SuOut1
(
    IOobject
    (
        "SuOut1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, Zero)
);

volScalarField SuOut2
(
    IOobject
    (
        "SuOut2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, Zero)
);

// Diffusion coefficient of acid in fluid phase
const dimensionedScalar diffInLiquid1
(
    "diffInLiquid",
    dimLength*dimLength/dimTime,
    solid1Inputs
);

const dimensionedScalar diffInLiquid2
(
    "diffInLiquid",
    dimLength*dimLength/dimTime,
    solid2Inputs
);

Info<< "Reading field diffusion coefficient for C1 field\n" << endl;
volScalarField diffCo1
(
    IOobject
    (
        "diffCo1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field diffusion coefficient for C2 field\n" << endl;
volScalarField diffCo2
(
    IOobject
    (
        "diffCo2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field for nucleation\n" << endl;
volScalarField nuSite
(
    IOobject
    (
        "nuSite",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volScalarField nuRate
(
    IOobject
    (
        "nuRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless/dimArea/dimTime, Zero)
);

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

dimensionedScalar D23("D23", dimViscosity, mixture);

// immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);
// 
// volScalarField& alpha1(mixture.alpha1());
// volScalarField& alpha2(mixture.alpha2());
// 
// const dimensionedScalar& rho1 = mixture.rho1();
// const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    // alpha1*rho1 + alpha2*rho2 + alpha3*rho3
    mesh,
    rho1
);
rho.oldTime();


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

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

dimensionedScalar liquidMass_(0.0);
dimensionedScalar liquidMassOld_(0.0);
dimensionedScalar solidMass_(0.0);
dimensionedScalar dmdtTot_(0.0);
dimensionedScalar deltaT_(0.0);

vectorList nuSiteList;
dimensionedScalar nuCount_(0.0);
dimensionedScalar nuTotal_(0.0);

scalar timeMarker = 0.0;

volScalarField Cmask
(
    IOobject
    (
        "Cmask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*0.0
);

phaseChangeReaction phaseChange1 = phaseChangeReaction(solid1Inputs, mesh, C1, alpha1, alpha2, Cmask);
phaseChangeReaction phaseChange2 = phaseChangeReaction(solid2Inputs, mesh, C2, alpha1, alpha3, Cmask);

//- diffusion coefficient flag and inputs
bool ionEffect = solid1Inputs.lookupOrDefault<bool>("ionEffect", false);
scalar ionEffectA = solid1Inputs.lookupOrDefault<scalar>("ionEffectA", 0.0);
scalar ionEffectB = solid1Inputs.lookupOrDefault<scalar>("ionEffectB", 0.0);
scalar conSat = solid1Inputs.lookupOrDefault<scalar>("Csat", 1.0);

//- initial nucleation vof
scalar nuIntVof = solid1Inputs.lookupOrDefault<scalar>("nucleationVof", 1e-6);
scalar nuIntVofCut = solid1Inputs.lookupOrDefault<scalar>("nucleationVofCut", 0.5);

//- criterion for removing solid phase from simulation domain
//scalar fluidOnly_ = solid1Inputs.lookupOrDefault<scalar>("fluidOnly", 0.99);

#include "createMRF.H"
#include "createFvOptions.H"
