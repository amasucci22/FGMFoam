/*-----------------------------------------------
Author:  Antonio Masucci
Created: 01/05/2026
Update:  gg/mm/yyyy

Progress Variable equation 
--------------------------------------------------*/


if (reacting)
{
	omegaYc = thermo.sourceYc();
}

else
{
	dimensionedScalar zero("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0);
	omegaYc = zero;
}

// Solve the Progress Variable equation

if (PDModel == "mixAvg")
{
	fvScalarMatrix YcEqn
    (
	  fvm::ddt(rho, Yc)
	  + fvm::div(phi, Yc)
	  - fvm::laplacian(thermo.dCC()+turbulence->mut()/Sct, Yc)
	  - fvc::laplacian(thermo.dCZ(), Z)	  
	  == 
	  omegaYc 
    );
   
    YcEqn.relax();
    YcEqn.solve();
    Yc = max(min( Yc, 1.0 ), 0.0 );
}

else if (PDModel == "Lewis")
{
	fvScalarMatrix YcEqn
    (
	  fvm::ddt(rho, Yc)
	  + fvm::div(phi, Yc)
	  - fvm::laplacian(turbulence->muEff()/Sct, Yc)
	  - fvc::laplacian(thermo.dCZ(), Z)	  
	  == 
	  omegaYc 
    );
   
    YcEqn.relax();
    YcEqn.solve();
    Yc = max(min( Yc, 1.0 ), 0.0 );
}
   
    
