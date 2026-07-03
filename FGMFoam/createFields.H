#include "createRDeltaT.H"
#include "readCombustionProperties.H"

Info<< "Reading thermophysical properties\n" << endl;


volScalarField ZV
(
    IOobject
    (
        "ZV",
        runTime.timeName(),
        mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField YcV
(
    IOobject
    (
        "YcV",
        runTime.timeName(),
        mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Yc
(
    IOobject
    (
        "Yc",
        runTime.timeName(),
        mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField YNO
(
    IOobject
    (
        "YNO",
        runTime.timeName(),
        mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField tangentialStrainRate
(
	IOobject
	(
        "tangentialStrainRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("tangentialStrainRate", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);
	
volScalarField stretch
(
	IOobject
	(
        "stretch",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("stretch", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);
	
volScalarField curvature
(
	IOobject
	(
        "curvature",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("curvature", dimensionSet(0,-1,0,0,0,0,0), 0.0)
);
	
volScalarField displacementSpeed
(
	IOobject
	(
        "displacementSpeed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("displacementSpeed", dimensionSet(0,1,-1,0,0,0,0), 0.0)
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

Info<< "Reading thermophysical properties\n" << endl;
// Initialize fgmTable and scalingYc table
lookupFGM FGMTable(mesh);
scaleYc scaleYcTable(mesh);

// Initializing fgmThermo
FGMThermo thermo(FGMModel,mesh,FGMTable,scaleYcTable,Yc,Z,YcV,ZV);
thermo.correct();

volScalarField& p = thermo.p();
//const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


volScalarField omegaYc
(
    IOobject
    (
        "omegaYc",
        runTime.timeName(),
        mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    thermo.sourceYc()
);

dimensionedScalar totalMass0=fvc::domainIntegrate(rho);
dimensionedScalar initialMass=fvc::domainIntegrate(rho);
p=totalMass0/fvc::domainIntegrate(psi);


#include "compressibleCreatePhi.H"


pressureReference pressureReference(p, pimple.dict(), thermo.incompressible());

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

volScalarField epsSGS
(
    IOobject
    (
        "epsSGS",
        runTime.timeName(),
        mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsSGS",dimensionSet(0,0,-1,0,0,0,0),1.e-6)
);

volScalarField kSGS
(
    IOobject
    (
        "kSGS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("kSGS",dimensionSet(0,2,-2,0,0,0,0),1e-6)
);

volScalarField deltaLES
(
    IOobject
    (
        "deltaLES",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedScalar("deltaLES",dimensionSet(0,1,0,0,0,0,0),1.00)
);

Info<< "Creating field beta_c\n" << endl;
volScalarField betac
(
	IOobject
		(
			"betaC",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		betacVal
);

volScalarField cellVolumes
(
    IOobject
    (
        "cellVolumes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVolume, 0.0)
);
cellVolumes.ref() = mesh.V();

if (!betacDynamic) betac = betacVal;

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
