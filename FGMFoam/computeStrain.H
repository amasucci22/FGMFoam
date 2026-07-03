//Calculate curvature, tangential strain rate and total stretch
volScalarField C = thermo.C();

// Compute the gradient of c (gradC is volVectorField)
volVectorField gradC = fvc::grad(C);

dimensionedScalar smallGradc("smallGradc", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1e-10);

// Normalize the gradient to get a unit vector (normgradc)
volVectorField normgradc = gradC / (mag(gradC) + smallGradc);

// Compute the curvature as the divergence of the normalized gradient (normgradc)
curvature = fvc::div(normgradc);

// Step 1: Compute s_n = divergence of U
volScalarField sn = fvc::div(U); // The divergence of the velocity field

// Step 2: Compute s_t = n dot grad(U) dot n (the tangential strain rate)
volScalarField st = normgradc & fvc::grad(U) & normgradc; // Compute n dot grad(U) dot n for each cell (scalar field)

// Step 3: Compute the tangential strain rate as sn - st
tangentialStrainRate = sn - st;

// Compute the displacement speed as the normalized material derivative of the normalized progress variable
displacementSpeed = 1/(mag(gradC) + smallGradc)*(fvc::ddt(C) + (U & gradC));

Info << "Computing Stretch" << endl;
// Compute the stretch
stretch = (curvature * displacementSpeed) + tangentialStrainRate;

/*volScalarField specificFlameSurface = sqrt(sqr(gradC.component(0)) + sqr(gradC.component(1)) + sqr(gradC.component(2)));

//scalar flameSurface = gSum(specificFlameSurface * cellVolumes);
// Force evaluation of the multiplication by creating an intermediate field
volScalarField integrand = specificFlameSurface * cellVolumes;

volScalarField deltaDiracPV = C - 0.8;

volScalarField deltaDiracPV_logical = C*0.0;


forAll(deltaDiracPV, cellI)
{
    if (mag(deltaDiracPV[cellI] < 0.01))
    {
        deltaDiracPV_logical[cellI] = 1;
    }
    else
    {
	deltaDiracPV_logical[cellI] = 0;
    }
}

// Force evaluation of the multiplication by creating an intermediate field
volScalarField integrandVervisch
(
    "integrandVervisch",
    deltaDiracPV_logical * mag(gradC) * cellVolumes //    specificFlameSurface * cellVolumes
);

// Now use gSum on the evaluated field
scalar flameSurfaceVervisch = gSum(integrandVervisch);

//scalar flameSurfaceVervisch = gSum(deltaDiracPV_logical & mag(gradC) & cellVolumes);

scalar consumptionSpeed = 0;

// Force evaluation of the multiplication by creating an intermediate field
volScalarField integrandSpeed
(
    "integrandSpeed",
    thermo.sourceH2() * cellVolumes //deltaDiracPV_logical * mag(gradC) * cellVolumes //    specificFlameSurface * cellVolumes
);

label inletPatchID = mesh.boundaryMesh().findPatchID("INLET");

const fvPatchScalarField& rhoB = rho.boundaryField()[inletPatchID];
const fvPatchScalarField& YFuelB = Z.boundaryField()[inletPatchID];

const scalarField& magSf = mesh.magSf().boundaryField()[inletPatchID];
scalar rhoU = gSum(rhoB * magSf) / gSum(magSf);
scalar YFuelU = gSum(YFuelB * magSf) / gSum(magSf);

consumptionSpeed = - gSum(integrandSpeed)/(flameSurfaceVervisch*YFuelU*rhoU);

volScalarField consumptionSpeedField
(
    IOobject
    (
        "consumptionSpeed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("consumptionSpeed", dimensionSet(0,1,-1,0,0,0,0), consumptionSpeed)
);

Info << "ConsumptionSpeed  = " << consumptionSpeed << " m/s" << endl;*/



