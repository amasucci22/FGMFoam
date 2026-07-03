{
	volScalarField& he = thermo.he();
	
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(phi, he)
      - fvm::laplacian(turbulence->muEff(), he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

}
