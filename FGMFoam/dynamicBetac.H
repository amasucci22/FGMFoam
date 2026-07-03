//Dynamic procedures must be called one time per timestep
//The following routines need to be modified with the tmp<...> format and object.clear() //Done
//Fix volScalarFields: Must be defined as IOobjects //Done

if (betacDynamic)
{  

   Info<<"Performing betac dynamic computation"<<endl;
   
   volScalarField C = thermo.C();
   volScalarField CV = thermo.CV();

   if (relTime==1)
   {
      thermo.correct();
      omegaYc = thermo.sourceYc();
   }
   
   kSGS = turbulence->k();
   dimensionedScalar ScLam("ScLam",dimless,0.7);
   dimensionedScalar bMin = 0.01;
   dimensionedScalar bMax = 100.0*betacVal;
   deltaLES = thermo.delta();

   Info<<"min/max betac: "<<min(betac).value()<<"/"<<max(betac).value()<<endl;

   volScalarField rhoFilt = thermo.filter(rho);
   volScalarField CFilt = thermo.filter(rho*C)/rhoFilt;

   volScalarField uDeltaFilt = thermo.filter(rho*sqrt(2.0*kSGS/3.0))/rhoFilt;
   volScalarField aij = thermo.C2(1.5*sqr(uDeltaFilt),2.0*T,2.0*deltaLES)*(thermo.filter(rho*CV));
   volScalarField bij = thermo.filter(rho*thermo.C2(kSGS,T,deltaLES)*CV);
   volScalarField Lij = thermo.filter(thermo.mu()/ScLam*magSqr(fvc::grad(C))) 
						- (thermo.filter(thermo.mu()/ScLam))*magSqr(fvc::grad(CFilt));

   const dimensionedScalar smallDim("smallDim",dimensionSet(2,-6,-2,0,0,0,0),1.0E-06);
   const dimensionedScalar smallDim1("smallDim1",dimensionSet(1,-3,-1,0,0,0,0),1.0E-06);
   aij.max(smallDim1*0.0);
   bij.max(smallDim1*0.0);

   
   scalarField& betaInt = betac.ref();
   const scalarField& aijInt = aij.internalField();
   const scalarField& bijInt = bij.ref();
   const scalarField& LijInt = Lij.internalField();
   const scalarField& CInt  = C.internalField();
   const scalarField& omInt  = omegaYc.internalField();
   const scalarField& rhoInt  = rho.internalField();

   scalar maxVal = 0.0;
   scalar minVal = 1.e6;

   // Updating betac values on the flame
   forAll (C.internalField(),celli)
   { 
      scalar& betai = betaInt[celli];
      const scalar& ai = aijInt[celli];
      const scalar& bi = bijInt[celli];
      const scalar& Li = LijInt[celli];
      const scalar& Ci = CInt[celli];
      const scalar& omegai = omInt[celli];
      const scalar& rhoi = rhoInt[celli];
      
      const scalar denom = Li; //ai-bi; //Li
 
      if ((omegai/rhoi>0.1) & (Ci>=0.05) & (Ci<=0.95) & (std::abs(denom)>1.e-6))
      {
             betai =  (ai-bi)/denom;
      }     
      
	  else //Need to avoid to put zero-dissipation when the SGS variance is small but non-zero
      {
         betai = 0.0; 
      }

              
      maxVal = max(maxVal,betai);
      minVal = min(minVal,betai);
   }
       
   reduce(minVal,minOp<scalar>());
   reduce(maxVal,maxOp<scalar>());

   Info<<"min/max betac: "<<min(betac).value()<<"/"<<max(betac).value()<<endl;
   Info<<"min/max betac Val: "<<minVal<<"/"<<maxVal<<endl;
     

   forAll(C.boundaryField(), patchi)
   {

       const fvPatchScalarField& CBound = C.boundaryField()[patchi];
       fvPatchScalarField& betaBound = betac.boundaryFieldRef()[patchi];

       forAll(CBound, facei)
       {
           label celli = mesh.boundary()[patchi].faceCells()[facei];
           betaBound[facei] = betaInt[celli];
       }
   }

   betac = thermo.filter(betac);
   betac.max(bMin);
   betac.min(bMax);

   rhoFilt.clear();
   CFilt.clear();
   uDeltaFilt.clear();
   aij.clear();
   bij.clear();
   Lij.clear();

}

