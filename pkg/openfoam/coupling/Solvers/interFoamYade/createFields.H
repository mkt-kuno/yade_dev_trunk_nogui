

#include "createRDeltaT.H"

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


Info<< "Reading transportProperties\n" << endl;

immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// const dimensionedScalar& nu1 = mixture.calcNu();
// const dimensionedScalar& nu1 = mixture.calcNu();


IOdictionary transportProperties
{
  IOobject
  (
   "transportProperties",
   runTime.constant(),
   mesh,
   IOobject::MUST_READ_IF_MODIFIED,
   IOobject::NO_WRITE )
};


//- YADE particle density
dimensionedScalar partDensity (
    "partDensity",
     dimDensity,
     transportProperties.lookup("partDensity")
);

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
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();



volScalarField mu = mixture.mu();

volScalarField nu = mu/rho;


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

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
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

Info << "Creating field alphac\n" << endl;
// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        "alphac",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);


volVectorField divT(
      IOobject(
      "divT",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
      fvc::laplacian(U)
);


volTensorField vGrad(
    IOobject(
      "vGrad",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      fvc::grad(U)
    );


volVectorField gradP(
      IOobject(
      "gradP",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE),
      fvc::grad(p)
);


volScalarField uSourceDrag
(
    IOobject
    (
      "uSourceDrag",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
    //   IOobject::NO_WRITE
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dragCoeff", dimensionSet(0,0,-1,0,0,0,0), 0 )
);

volScalarField uCoeff
(
    IOobject
    (
        "uCoeffDrag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        // IOobject::NO_WRITE
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dragCoeffFluid", dimensionSet(0,0,-1,0,0,0,0), 0 
    )
);


volVectorField uInterp(
    IOobject(
      "uInterpField",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
    mesh,
    dimensionedVector("interpVel", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) )
);


volVectorField uSource(
    IOobject(
      "uSource",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE),
    //   IOobject::NO_WRITE),
    mesh,
    dimensionedVector("volSource", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0) )
);

volVectorField uParticle(
    IOobject(
      "uParticle",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
    //   IOobject::NO_WRITE),
      IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("uParticle", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) )
);


volVectorField ddtU_f (
    IOobject (
      "ddtU_f",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE),
    //   IOobject::NO_WRITE),
      fvc::ddt(U)+fvc::div(phi, U)
);



surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));

surfaceScalarField alphaPhic("alphaPhic", alphacf*phi);


#include "createMRF.H"
#include "createFvOptions.H"
