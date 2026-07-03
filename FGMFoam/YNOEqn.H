/*-----------------------------------------------
Author:  Antonio Masucci
Created: 01/05/2026
Update:  gg/mm/yyyy

NO Transport equation 
--------------------------------------------------*/


if (PDModel == "mixAvg")
{
	fvScalarMatrix YNOEqn
	(
	  fvm::ddt(rho, YNO)
	  + fvm::div(phi, YNO)
	  - fvm::laplacian(thermo.rho()*thermo.DNO(), YNO)  
	  == 
	  thermo.sourceNO() 
	);

	YNOEqn.relax();   
	YNOEqn.solve();
	YNO = max(min( YNO, 1.0 ), 0.0 );
}

if (PDModel == "Lewis")
{
	fvScalarMatrix YNOEqn
	(
	  fvm::ddt(rho, YNO)
	  + fvm::div(phi, YNO)
	  - fvm::laplacian(turbulence->muEff()/Sct, YNO)  
	  == 
	  thermo.sourceNO() 
	);

	YNOEqn.relax();   
	YNOEqn.solve();
	YNO = max(min( YNO, 1.0 ), 0.0 );
}

    
