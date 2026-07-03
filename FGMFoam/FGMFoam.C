/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    FGMFoam

Description
    Transient solver for turbulent flow of Low Mach reacting fluids with
	FGM model and Beta-FDF approach.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "lookupFGM.H"
#include "scaleYc.H"
#include "FGMThermo.H"
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createRhoUfIfPresent.H"
	
	turbulence->validate();


    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
	int relTime = 0;

    while (runTime.run())
    {
		#include "readDyMControls.H"
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

		
		// Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
		autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }
		
		runTime++;
		relTime++;
		
		Info<< nl << "Time = " << runTime.timeName() << endl;
		Info<< "relTime = " << relTime << nl << endl;
		
		//#include "dynamicBetac.H"

        while (pimple.loop())
        {
			
			if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", rho*U);
                }

                fvModels.preUpdateMesh();

                // Do any mesh changes
                mesh.update();
			}
			
			if (
             !mesh.steady()
             && !pimple.simpleRho()
             && pimple.firstPimpleIter()
            )
            {													  
                #include "rhoEqn.H"
            }

            fvModels.correct();
			
			kSGS = turbulence->k();
			#include "UEqn.H"
			#include "ZEqn.H"
			#include "ZVEqn.H"
			#include "YcEqn.H"
			#include "YcVEqn.H"
			#include "YNOEqn.H"
			thermo.correct();


            // --- Pressure corrector loop
            while (pimple.correct())
            {
               #include "pEqn.H"              
            }
			
			if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }

        }

        if (!mesh.steady())
        {
            rho = thermo.rho();
        }
		
		//#include "computeStrain.H"

		Info<< "Yc min/max   = " << min(Yc).value() << "/" << max(Yc).value() << endl;
		Info<< "Z min/max    = " << min(Z).value() << "/" << max(Z).value() << endl;
		Info<< "CV min/max   = " << min(thermo.CV()).value() << "/" << max(thermo.CV()).value() << endl;
		Info<< "ZV min/max   = " << min(ZV).value() << "/" << max(ZV).value() << endl;
		Info<< "SYc min/max  = " << min(thermo.sourceYc()).value() << "/" << max(thermo.sourceYc()).value() << endl;
		Info<< "SYcV min/max = " << min(thermo.sourceYcV()).value() << "/" << max(thermo.sourceYcV()).value() << endl;
		Info<< "T min/max    = " << min(thermo.T()).value() << "/" << max(thermo.T()).value() << endl;
		Info<< "rho min/max  = " << min(rho).value() << "/" << max(rho).value() << endl;
		Info<< "p min/max    = " << min(p).value() << "/ " << max(p).value() << endl;
		Info<< "WMix min/max = " << min(thermo.W()).value() << "/ " << max(thermo.W()).value() << endl;
		Info<< "U min/max    = " << min(U).value() << "/" << max(U).value() << endl;
		
        runTime.write();

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
