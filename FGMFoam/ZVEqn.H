/*-----------------------------------------------
Author:  Antonio Masucci
Created: 01/05/2026
Update:  gg/mm/yyyy

Mixture Fraction Variance equation 
--------------------------------------------------*/

//volScalarField zvSGSdiss = thermo.zvSGSdiss(turbulence->mut(),thermo.delta())*ZV/(0.5*Sct);
//epsSGS = thermo.C2(kSGS,T,thermo.delta())*YcV;
volScalarField S = symm(fvc::grad(U)) && dev(symm(fvc::grad(U)));

const dimensionedScalar smallDim3("smallDim3",dimensionSet(0,2,-3,0,0,0,0),1.0E-20);
volScalarField tauDelta = pow((magSqr(1.8)* turbulence->nut() * S + smallDim3)/magSqr(thermo.delta()),-1.0/3.0);

volScalarField tau_turb = pow
(
    (magSqr(1.8)*turbulence->nut()*S + smallDim3)
    / magSqr(thermo.delta()),
    -1.0/3.0
);

volScalarField tau_lam = sqr(thermo.delta()) / (thermo.dCC()/thermo.rho());

// Choose between tauTurb and tauLam do dissipate YcV
tauDelta = min(tau_turb, tau_lam);

// Solve the Mixture Fraction Variance equation
if (PDModel == "mixAvg")
{
	fvScalarMatrix ZVEqn
    (
	  fvm::ddt(rho, ZV)
	  + fvm::div(phi, ZV)
	  - fvm::laplacian(thermo.dZZ()+turbulence->mut()/Sct, ZV)
	  - 2.0*turbulence->mut()/Sct*magSqr(fvc::grad(Z))
	  + 2.0*rho*ZV/tauDelta
	  //+ fvm::Sp(2.0*thermo.zvSGSdiss(turbulence->mut(),thermo.delta())/(Sct*0.5),ZV)	  
    );
   
	ZVEqn.relax();
    ZVEqn.solve();
    ZV = max(min( ZV, 0.25 ), 0.0 );
}

else if (PDModel == "Lewis")
{
	fvScalarMatrix ZVEqn
    (
	  fvm::ddt(rho, ZV)
	  + fvm::div(phi, ZV)
	  - fvm::laplacian(turbulence->muEff()/Sct, ZV)
	  - 2.0*turbulence->mut()/Sct*magSqr(fvc::grad(Z))
	  + 2.0*rho*ZV/tauDelta
	  //+ fvm::Sp(2.0*thermo.zvSGSdiss(turbulence->mut(),thermo.delta())/(Sct*0.5),ZV)	  
    );
   
	ZVEqn.relax();
    ZVEqn.solve();
    ZV = max(min( ZV, 0.25 ), 0.0 );
}



    
