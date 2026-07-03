/*-----------------------------------------------
Author:  Antonio Masucci
Created: 01/05/2026
Update:  gg/mm/yyyy

Progress Variable Variance equation 
--------------------------------------------------*/


if (PDModel == "mixAvg")
{
	fvScalarMatrix YcVEqn
	(
	  fvm::ddt(rho, YcV)
	  + fvm::div(phi, YcV)
	  - fvm::laplacian(thermo.dCC()+turbulence->mut()/Sct, YcV)
	  - 2.0*turbulence->mut()/Sct*magSqr(fvc::grad(Yc))
	  + 2.0*rho*YcV/tauDelta
	  //+ fvm::Sp(2.0*rho*thermo.C2(kSGS,T,thermo.delta())/betac,YcV)	 	  
	  == 
	  2.0*(thermo.sourceYcV()-Yc*omegaYc) 
    );
   
	YcVEqn.relax();
    YcVEqn.solve();
    YcV = max(min( YcV, 0.25 ), 0.0 );
}

else if (PDModel == "Lewis")
{
	fvScalarMatrix YcVEqn
	(
	  fvm::ddt(rho, YcV)
	  + fvm::div(phi, YcV)
	  - fvm::laplacian(turbulence->muEff()/Sct, YcV)
	  - 2.0*turbulence->mut()/Sct*magSqr(fvc::grad(Yc))
	  + 24.0*rho*YcV/tauDelta
	  //+ fvm::Sp(2.0*rho*thermo.C2(kSGS,T,thermo.delta())/betac,YcV)	 	  
	  == 
	  2.0*(thermo.sourceYcV()-Yc*omegaYc) 
    );
   
	YcVEqn.relax();
    YcVEqn.solve();
    YcV = max(min( YcV, 0.25 ), 0.0 );
}


    
