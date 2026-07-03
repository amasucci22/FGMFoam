/*-----------------------------------------------
Author:  Antonio Masucci
Created: 01/05/2026
Update:  gg/mm/yyyy

Mixture Fraction equation 
--------------------------------------------------*/
const scalar Sct = 0.7;

if (PDModel == "mixAvg")
{
	fvScalarMatrix ZEqn
	(
	  fvm::ddt(rho, Z)
	  + fvm::div(phi, Z)
	  - fvm::laplacian(thermo.dZZ()+turbulence->mut()/Sct, Z)
	  - fvc::laplacian(thermo.dZC(), Yc)  
    );
	
	ZEqn.relax();
    ZEqn.solve();
    Z = max(min( Z, 1.0 ), 0.0 );
}

else if (PDModel == "Lewis")
{
	fvScalarMatrix ZEqn
	(
	  fvm::ddt(rho, Z)
	  + fvm::div(phi, Z)
	  - fvm::laplacian(turbulence->muEff()/Sct, Z)
	  - fvc::laplacian(thermo.dZC(), Yc)  
    );
	
	ZEqn.relax();
    ZEqn.solve();
    Z = max(min( Z, 1.0 ), 0.0 );
}


    
