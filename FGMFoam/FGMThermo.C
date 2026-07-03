/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "FGMThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "fvm.H"
#include "lookupFGM.H"
#include "scaleYc.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

namespace Foam
{
	
void FGMThermo::heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hBf[patchi]).gradient()
                = hBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hBf[patchi]).refGrad()
                = hBf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FGMThermo::FGMThermo
(
	const word& model,
    const fvMesh& mesh,
    const lookupFGM& lookupFGM,
	const scaleYc& scaleYc,
    volScalarField& Yc,
	volScalarField& Z,
	volScalarField& YcV,
	volScalarField& ZV,
	const word& phaseName
)
:
    psiThermo::composite(mesh,phaseName),

	FGMTable_(lookupFGM),
    scaleYcTable_(scaleYc),
	
	Yc_(Yc),
    Z_(Z),
	YcV_(YcV),
    ZV_(ZV),
	
    he_
    (
        IOobject
        (
            "he",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        this->heBoundaryTypes(),  
        this->heBoundaryBaseTypes()
    ),

    Cp_
    (
    IOobject
    (
        "Cp",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(FGMTable_.lookup("CpAir"))
    ),

    Cv_
   (
    IOobject
    (
        "Cv",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(FGMTable_.lookup("CpAir"))
    ),
	
	dH_
   (
    IOobject
    (
        "dH",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(FGMTable_.lookup("dHAir"))
    ),
	
	Temp_
   (
       IOobject
       (
           "T",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionedScalar(FGMTable_.lookup("T0"))
    ),
	
	YcMin_
   (
       IOobject
       (
           "YcMin",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionSet(0,0,0,0,0,0,0)
    ),
	
	YcMax_
   (
       IOobject
       (
           "YcMax",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionedScalar("YcMax",dimensionSet(0,0,0,0,0,0,0),1e-6)
    ),
	
	ZScal_
   (
       IOobject
       (
           "ZScal",
           mesh.time().timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionSet(0,0,0,0,0,0,0)
    ),

	C_
   (
    IOobject
    (
        "C",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("C", dimless, 0)
   ),
  
  sourceYc_
  (
    IOobject
    (
        "sourceYc",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1,-3,-1,0,0,0,0)
   ),
   
   sourceYcV_
  (
    IOobject
    (
        "sourceYcV",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1,-3,-1,0,0,0,0)
   ),
   
   sourceH2_
  (
    IOobject
    (
        "sourceH2",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1,-3,-1,0,0,0,0)
   ),
   
   sourceNO_
  (
    IOobject
    (
        "sourceNO",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1,-3,-1,0,0,0,0)
   ),
   
  dCC_
  (
    IOobject
    (
        "dCC",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume*dimLength*dimLength/dimTime,1E-5)
   ),
   
  dCZ_
  (
    IOobject
    (
        "dCZ",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume*dimLength*dimLength/dimTime,1E-5)
   ),
   
   dZC_
  (
    IOobject
    (
        "dZC",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume*dimLength*dimLength/dimTime,1E-5)
   ),
   
   dZZ_
  (
    IOobject
    (
        "dZZ",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume*dimLength*dimLength/dimTime,1E-5)
   ),
   
   k_
  (
    IOobject
    (
        "k",
        mesh.time().timeName(),
        mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimTime/dimLength/dimTemperature,1E-5)
   ),
   
   D_
   (
    IOobject
    (
        "D",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume*dimLength*dimLength/dimTime,1E-5)
    ),
	
	 DNO_
   (
    IOobject
    (
        "DNO",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimLength/dimTime,1E-5)
    ),
	W_
   (
    IOobject
    (
        "W",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(FGMTable_.lookup("WAir"))
    ),
	
	H_
   (
    IOobject
    (
        "H",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	H2_
   (
    IOobject
    (
        "H2",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	H2O_
   (
    IOobject
    (
        "H2O",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	HO2_
   (
    IOobject
    (
        "HO2",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	H2O2_
   (
    IOobject
    (
        "H2O2",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	O_
   (
    IOobject
    (
        "O",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	O2_
   (
    IOobject
    (
        "O2",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	OH_
   (
    IOobject
    (
        "OH",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	N2_
   (
    IOobject
    (
        "N2",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	NO_
   (
    IOobject
    (
        "NO",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	HRR_
   (
    IOobject
    (
        "HRR",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1,-1,-3,0,0,0,0)
    ),
	
	sL_
   (
    IOobject
    (
        "sL",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,1,-1,0,0,0,0)
    ),
	
	lF_
   (
    IOobject
    (
        "lF",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,1,0,0,0,0,0)
    ),
	
	tau_
   (
    IOobject
    (
        "tau",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0)
    ),
	
	CV_
   (
    IOobject
    (
        "CVar",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CVar",dimensionSet(0,0,0,0,0,0,0),0.0)
    )
	
{
	correct();
	
	const dimensionedScalar T0(FGMTable_.lookup("T0"));
	
	Info<< "T0 = " << T0.value() << endl;
	
    he_ == he(this->p_,this->T_);
	
	Info<< "T0 = " << T0.value() << endl;
	
    heBoundaryCorrection(he_);
	
	this->psi_.oldTime();
	
	
}




// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FGMThermo::~FGMThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FGMThermo::correct()
{
  scalarField& TCells 			= Temp_.primitiveFieldRef();
  scalarField& CCells 			= C_.primitiveFieldRef();
  scalarField& YcMinCells 		= YcMin_.primitiveFieldRef();
  scalarField& YcMaxCells 		= YcMax_.primitiveFieldRef();
  scalarField& sLCells 			= sL_.primitiveFieldRef();
  scalarField& lFCells 			= lF_.primitiveFieldRef();
  scalarField& tauCells 		= tau_.primitiveFieldRef();
  scalarField& dHCells 			= dH_.primitiveFieldRef();
  scalarField& sourceYcCells 	= sourceYc_.primitiveFieldRef();
  scalarField& sourceYcVCells 	= sourceYcV_.primitiveFieldRef();
  scalarField& sourceH2Cells 	= sourceH2_.primitiveFieldRef();
  scalarField& sourceNOCells    = sourceNO_.primitiveFieldRef();

  scalarField& dCCCells 		= dCC_.primitiveFieldRef();
  scalarField& dCZCells 		= dCZ_.primitiveFieldRef();
  scalarField& dZCCells 		= dZC_.primitiveFieldRef();
  scalarField& dZZCells 		= dZZ_.primitiveFieldRef();
  scalarField& kCells 		= k_.primitiveFieldRef();
  scalarField& CpCells 			= Cp_.primitiveFieldRef();
  scalarField& WCells 			= W_.primitiveFieldRef();
  scalarField& DCells 			= D_.primitiveFieldRef();
  scalarField& DNOCells         = DNO_.primitiveFieldRef();
									   
  scalarField& HCells 			= H_.primitiveFieldRef();
  scalarField& H2Cells 			= H2_.primitiveFieldRef();
  scalarField& H2OCells 		= H2O_.primitiveFieldRef();
  scalarField& HO2Cells 		= HO2_.primitiveFieldRef();
  scalarField& H2O2Cells 		= H2O2_.primitiveFieldRef();
  scalarField& OHCells 			= OH_.primitiveFieldRef();
  scalarField& OCells 			= O_.primitiveFieldRef();
  scalarField& O2Cells 			= O2_.primitiveFieldRef();
  scalarField& N2Cells 			= N2_.primitiveFieldRef();
  scalarField& NOCells          = NO_.primitiveFieldRef();

  scalarField& HRRCells 		= HRR_.primitiveFieldRef();
  scalarField& CVCells 			= CV_.primitiveFieldRef();
  
  const scalarField& YcCells 	= Yc_.internalField();
  const scalarField& ZCells 	= Z_.internalField();
  const scalarField& YcVCells 	= YcV_.internalField();
  const scalarField& ZVCells 	= ZV_.internalField();
  
  const int nPointsC 			= FGMTable_.nC_;
  const int nPointsCv 			= FGMTable_.nCv_;
  const int nPointsZv 			= FGMTable_.nZv_;
  
  const scalarList& ZTab 		= FGMTable_.Z_table;
  const scalarList& CTab 		= FGMTable_.C_table;
  const scalarList& ZVTab 		= FGMTable_.ZV_table;
  const scalarList& CVTab 		= FGMTable_.CV_table;
  const scalarList& sLTab 		= FGMTable_.sL_table;
  const scalarList& lFTab 		= FGMTable_.lF_table;
  const scalarList& tauTab 		= FGMTable_.tau_table;
  const scalarList& sourceYcTab = FGMTable_.sourceYc_table;
  const scalarList& sourceYcVTab= FGMTable_.sourceYcV_table;
  const scalarList& sourceH2Tab = FGMTable_.sourceH2_table;
  const scalarList& sourceNOTab = FGMTable_.sourceNO_table;
  const scalarList& TTab 		= FGMTable_.Temp_table;
  const scalarList& kTab 	= FGMTable_.k_table;
  const scalarList& dHTab 		= FGMTable_.dH_table;
  const scalarList& CpTab 		= FGMTable_.Cp_table;
  const scalarList& DTab 		= FGMTable_.D_table;
  const scalarList& DNOTab 		= FGMTable_.DNO_table;
  const scalarList& WTab 		= FGMTable_.W_table;
  const scalarList& HTab 		= FGMTable_.H_table;
  const scalarList& H2Tab 		= FGMTable_.H2_table;
  const scalarList& H2OTab 		= FGMTable_.H2O_table;
  const scalarList& H2O2Tab 	= FGMTable_.H2O2_table;
  const scalarList& HO2Tab 		= FGMTable_.HO2_table;
  const scalarList& OTab 		= FGMTable_.O_table;
  const scalarList& O2Tab 		= FGMTable_.O2_table;
  const scalarList& OHTab 		= FGMTable_.OH_table;
  const scalarList& N2Tab 		= FGMTable_.N2_table;
  const scalarList& NOTab 		= FGMTable_.NO_table;
  const scalarList& HRRTab 		= FGMTable_.HRR_table;
  const scalarList& dCCTab 		= FGMTable_.dCC_table;
  const scalarList& dCZTab 		= FGMTable_.dCZ_table;
  const scalarList& dZCTab 		= FGMTable_.dZC_table;
  const scalarList& dZZTab 		= FGMTable_.dZZ_table;
  
  const scalarList& YcMinTab 	= scaleYcTable_.YcMin;
  const scalarList& YcMaxTab 	= scaleYcTable_.YcMax;
  const scalarList& ZScalTab 	= scaleYcTable_.ZScal;
  
  const double zLean 			= FGMTable_.zLean_.value();
  const double zRich 			= FGMTable_.zRich_.value();
  
  const dimensionedScalar T0(FGMTable_.lookup("T0"));
  const dimensionedScalar R0(FGMTable_.lookup("R0"));
  
  
  
  // Interpolate for internal field
  forAll(ZCells, celli)
	  {
		  YcMinCells[celli] = scaleYcTable_.interpolateValue1D(YcMinTab,ZCells[celli],ZScalTab);
		  YcMaxCells[celli] = scaleYcTable_.interpolateValue1D(YcMaxTab,ZCells[celli],ZScalTab);										   
		  CCells[celli] 	= (YcCells[celli]-YcMinCells[celli])/(YcMaxCells[celli]-YcMinCells[celli]+SMALL);
		  CCells[celli] 	= max(min( CCells[celli], 1.0 ), 0.0 );
		  CVCells[celli] 	= (YcVCells[celli])/(magSqr(YcMaxCells[celli])+1E-12);
		  CVCells[celli] 	= max(min( CVCells[celli], 0.25 ), 0.0 );
		  
		  const double& Zi  		= ZCells[celli];
		  const double& ZVi 		= ZVCells[celli];
		  const double& Ci  		= CCells[celli];
		  const double& CVi 		= CVCells[celli];
		  

		  if (Zi < zLean)
		  { 			 
			const double ZEff 		= ZCells[celli]/zLean;
			
			CCells[celli]     		= CCells[celli]*ZEff;
			CVCells[celli]    		= CVCells[celli]*ZEff;
			CCells[celli] 			= max(min( CCells[celli], 1.0 ), 0.0 );
			CVCells[celli] 			= max(min( CVCells[celli], 0.25 ), 0.0 );
			
			const double& Ci  		= CCells[celli];
			const double& CVi 		= CVCells[celli];
			
			sourceYcCells[celli] 	= FGMTable_.interpolateValue4D(sourceYcTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceYcVCells[celli] 	= FGMTable_.interpolateValue4D(sourceYcVTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceH2Cells[celli] 	= FGMTable_.interpolateValue4D(sourceH2Tab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceNOCells[celli] 	= FGMTable_.interpolateValue4D(sourceNOTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sLCells[celli] 			= FGMTable_.interpolateValue1D(sLTab,zLean,ZTab);
			lFCells[celli] 			= FGMTable_.interpolateValue1D(lFTab,zLean,ZTab);
			tauCells[celli] 		= FGMTable_.interpolateValue1D(tauTab,zLean,ZTab);
			TCells[celli] 			= FGMTable_.interpolateValue4D(TTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
		    CpCells[celli] 			= FGMTable_.interpolateValue4D(CpTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			WCells[celli] 			= FGMTable_.interpolateValue4D(WTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			DCells[celli] 			= FGMTable_.interpolateValue4D(DTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			DNOCells[celli] 		= FGMTable_.interpolateValue4D(DNOTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			kCells[celli]  	        = FGMTable_.interpolateValue4D(kTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dHCells[celli] 			= FGMTable_.interpolateValue4D(dHTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HCells[celli] 			= FGMTable_.interpolateValue4D(HTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2Cells[celli] 			= FGMTable_.interpolateValue4D(H2Tab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2OCells[celli] 		= FGMTable_.interpolateValue4D(H2OTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2O2Cells[celli] 		= FGMTable_.interpolateValue4D(H2O2Tab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HO2Cells[celli] 		= FGMTable_.interpolateValue4D(HO2Tab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			OCells[celli] 			= FGMTable_.interpolateValue4D(OTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			O2Cells[celli] 			= FGMTable_.interpolateValue4D(O2Tab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			OHCells[celli] 			= FGMTable_.interpolateValue4D(OHTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			N2Cells[celli] 			= FGMTable_.interpolateValue4D(N2Tab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			NOCells[celli] 			= FGMTable_.interpolateValue4D(NOTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dCCCells[celli] 		= FGMTable_.interpolateValue4D(dCCTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dCZCells[celli] 		= FGMTable_.interpolateValue4D(dCZTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dZCCells[celli] 		= FGMTable_.interpolateValue4D(dZCTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dZZCells[celli] 		= FGMTable_.interpolateValue4D(dZZTab,zLean,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			
			
			TCells[celli]     		= FGMTable_.T0_.value() + (TCells[celli] - FGMTable_.T0_.value())*ZEff;
			sourceYcCells[celli]    = 0.0;
			sourceYcVCells[celli]   = 0.0;
			sourceH2Cells[celli]    = 0.0;
			sourceNOCells[celli]    = sourceNOCells[celli]*ZEff;
			HRRCells[celli]         = HRRCells[celli]*ZEff;
			sLCells[celli]    		= sLCells[celli]*ZEff;
			lFCells[celli]    		= lFCells[celli]*ZEff;
			tauCells[celli]    		= tauCells[celli]*ZEff;
			CpCells[celli]      	= FGMTable_.CpAir_.value() + (CpCells[celli] - FGMTable_.CpAir_.value())*ZEff;
			WCells[celli]     		= FGMTable_.WAir_.value() + (WCells[celli] - FGMTable_.WAir_.value())*ZEff;
			DCells[celli]     		= FGMTable_.DAir_.value() + (DCells[celli] - FGMTable_.DAir_.value())*ZEff;
			DNOCells[celli]         = FGMTable_.DNOAir_.value() + (DNOCells[celli] - FGMTable_.DNOAir_.value())*ZEff;
			dHCells[celli]     		= FGMTable_.dHAir_.value() + (dHCells[celli] - FGMTable_.dHAir_.value())*ZEff;
			kCells[celli]  	        = FGMTable_.kAir_.value() + (kCells[celli] - FGMTable_.kAir_.value())*ZEff;
			HCells[celli]     		= FGMTable_.HAir_.value() + (HCells[celli] - FGMTable_.HAir_.value())*ZEff;
			H2Cells[celli]     		= FGMTable_.H2Air_.value() + (H2Cells[celli] - FGMTable_.H2Air_.value())*ZEff;
			H2OCells[celli]     	= FGMTable_.H2OAir_.value() + (H2OCells[celli] - FGMTable_.H2OAir_.value())*ZEff;
			H2O2Cells[celli]    	= FGMTable_.H2O2Air_.value() + (H2O2Cells[celli] - FGMTable_.H2O2Air_.value())*ZEff;
			HO2Cells[celli]     	= FGMTable_.HO2Air_.value() + (HO2Cells[celli] - FGMTable_.HO2Air_.value())*ZEff;
			OCells[celli]     		= FGMTable_.OAir_.value() + (OCells[celli] - FGMTable_.OAir_.value())*ZEff;
			O2Cells[celli]     		= FGMTable_.O2Air_.value() + (O2Cells[celli] - FGMTable_.O2Air_.value())*ZEff;
			OHCells[celli]     		= FGMTable_.OHAir_.value() + (OHCells[celli] - FGMTable_.OHAir_.value())*ZEff;
			N2Cells[celli]     		= FGMTable_.N2Air_.value() + (N2Cells[celli] - FGMTable_.N2Air_.value())*ZEff;
			NOCells[celli]     		= NOCells[celli]*ZEff;
			dCCCells[celli]     	= FGMTable_.dCCAir_.value() + (dCCCells[celli] - FGMTable_.dCCAir_.value())*ZEff;
			dCZCells[celli]     	= FGMTable_.dCZAir_.value() + (dCZCells[celli] - FGMTable_.dCZAir_.value())*ZEff;
			dZCCells[celli]     	= FGMTable_.dZCAir_.value() + (dZCCells[celli] - FGMTable_.dZCAir_.value())*ZEff;
			dZZCells[celli]     	= FGMTable_.dZZAir_.value() + (dZZCells[celli] - FGMTable_.dZZAir_.value())*ZEff; 
		  }
		  
		  else if (Zi > zRich)
		  { 		 
			const double ZEff 		= (Zi-zRich)/(1-zRich);
			
			CCells[celli]     		= CCells[celli] - CCells[celli]*ZEff;
			CVCells[celli]     		= CVCells[celli] - CVCells[celli]*ZEff;
			CCells[celli] 			= max(min( CCells[celli], 1.0 ), 0.0 );
			CVCells[celli] 			= max(min( CVCells[celli], 0.25 ), 0.0 );
			
			const double& Ci  		= CCells[celli];
		    const double& CVi 		= CVCells[celli];
			
			TCells[celli] 			= FGMTable_.interpolateValue4D(TTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceNOCells[celli] 	= FGMTable_.interpolateValue4D(sourceNOTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
		    CpCells[celli] 			= FGMTable_.interpolateValue4D(CpTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			WCells[celli] 			= FGMTable_.interpolateValue4D(WTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			DCells[celli] 			= FGMTable_.interpolateValue4D(DTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			DNOCells[celli] 		= FGMTable_.interpolateValue4D(DNOTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			kCells[celli] 		    = FGMTable_.interpolateValue4D(kTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dHCells[celli] 			= FGMTable_.interpolateValue4D(dHTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HCells[celli] 			= FGMTable_.interpolateValue4D(HTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2Cells[celli] 			= FGMTable_.interpolateValue4D(H2Tab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2OCells[celli] 		= FGMTable_.interpolateValue4D(H2OTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2O2Cells[celli] 		= FGMTable_.interpolateValue4D(H2O2Tab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HO2Cells[celli] 		= FGMTable_.interpolateValue4D(HO2Tab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			OCells[celli] 			= FGMTable_.interpolateValue4D(OTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			O2Cells[celli] 			= FGMTable_.interpolateValue4D(O2Tab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			OHCells[celli] 			= FGMTable_.interpolateValue4D(OHTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			N2Cells[celli] 			= FGMTable_.interpolateValue4D(N2Tab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			NOCells[celli] 			= FGMTable_.interpolateValue4D(NOTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dCCCells[celli] 		= FGMTable_.interpolateValue4D(dCCTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dCZCells[celli] 		= FGMTable_.interpolateValue4D(dCZTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dZCCells[celli] 		= FGMTable_.interpolateValue4D(dZCTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dZZCells[celli] 		= FGMTable_.interpolateValue4D(dZZTab,zRich,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			
			sourceYcCells[celli]    = 0.0;
			sourceYcVCells[celli]   = 0.0;
			sourceH2Cells[celli]    = 0.0;
			HRRCells[celli]     	= 0.0;
			sLCells[celli] 			= 0.0;
			lFCells[celli] 			= 0.0;
			tauCells[celli] 		= 0.0;
			
			TCells[celli]     		= TCells[celli] + (FGMTable_.T0_.value() - TCells[celli])*ZEff;
			sourceNOCells[celli]    = sourceNOCells[celli] - sourceNOCells[celli]*ZEff;
			CpCells[celli]     		= CpCells[celli] + (FGMTable_.CpFuel_.value() - CpCells[celli])*ZEff;
			WCells[celli]     		= WCells[celli] + (FGMTable_.WFuel_.value() - WCells[celli])*ZEff;
			DCells[celli]     		= DCells[celli] + (FGMTable_.DFuel_.value() - DCells[celli])*ZEff;
			DNOCells[celli]     	= DNOCells[celli] + (FGMTable_.DNOFuel_.value() - DNOCells[celli])*ZEff;
			dHCells[celli]     		= dHCells[celli] + (FGMTable_.dHFuel_.value() - dHCells[celli])*ZEff;
			kCells[celli]  	        = kCells[celli] + (FGMTable_.kFuel_.value() - kCells[celli])*ZEff;
			HCells[celli]     		= HCells[celli] + (FGMTable_.HFuel_.value() - HCells[celli])*ZEff;
			H2Cells[celli]     		= H2Cells[celli] + (FGMTable_.H2Fuel_.value() - H2Cells[celli])*ZEff;
			H2OCells[celli]     	= H2OCells[celli] + (FGMTable_.H2OFuel_.value() - H2OCells[celli])*ZEff;
			H2O2Cells[celli]    	= FGMTable_.H2O2Fuel_.value() + (H2O2Cells[celli] - FGMTable_.H2O2Fuel_.value())*ZEff;
			HO2Cells[celli]     	= HO2Cells[celli] + (FGMTable_.HO2Fuel_.value() - HO2Cells[celli])*ZEff;
			OCells[celli]     		= OCells[celli] + (FGMTable_.OFuel_.value() - OCells[celli])*ZEff;
			O2Cells[celli]     		= O2Cells[celli] + (FGMTable_.O2Fuel_.value() - O2Cells[celli])*ZEff;
			OHCells[celli]     		= OHCells[celli] + (FGMTable_.OHFuel_.value() - OHCells[celli])*ZEff;
			N2Cells[celli]     		= N2Cells[celli] + (FGMTable_.N2Fuel_.value() - N2Cells[celli])*ZEff;
			NOCells[celli]     		= NOCells[celli] - NOCells[celli]*ZEff;
			dCCCells[celli]     	= dCCCells[celli] + (FGMTable_.dCCFuel_.value() - dCCCells[celli])*ZEff;
			dCZCells[celli]     	= dCZCells[celli] + (FGMTable_.dCZFuel_.value() - dCZCells[celli])*ZEff;
			dZCCells[celli]     	= dZCCells[celli] + (FGMTable_.dZCFuel_.value() - dZCCells[celli])*ZEff;
			dZZCells[celli]     	= dZZCells[celli] + (FGMTable_.dZZFuel_.value() - dZZCells[celli])*ZEff;
		  }
				
		  else
		  {
			sourceYcCells[celli] 	= FGMTable_.interpolateValue4D(sourceYcTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceYcVCells[celli] 	= FGMTable_.interpolateValue4D(sourceYcVTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceH2Cells[celli] 	= FGMTable_.interpolateValue4D(sourceH2Tab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceNOCells[celli] 	= FGMTable_.interpolateValue4D(sourceNOTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			sourceYcCells[celli] 	= max(sourceYcCells[celli], 0.0 );
			sourceYcVCells[celli]	= max(sourceYcVCells[celli], 0.0 );
			  
			sLCells[celli] 			= FGMTable_.interpolateValue1D(sLTab,Zi,ZTab);
			lFCells[celli] 			= FGMTable_.interpolateValue1D(lFTab,Zi,ZTab);
			tauCells[celli] 		= FGMTable_.interpolateValue1D(tauTab,Zi,ZTab);
			  
			TCells[celli] 			= FGMTable_.interpolateValue4D(TTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			kCells[celli] 		    = FGMTable_.interpolateValue4D(kTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dHCells[celli] 			= FGMTable_.interpolateValue4D(dHTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			CpCells[celli] 			= FGMTable_.interpolateValue4D(CpTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			WCells[celli] 			= FGMTable_.interpolateValue4D(WTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			DCells[celli] 			= FGMTable_.interpolateValue4D(DTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			DNOCells[celli] 		= FGMTable_.interpolateValue4D(DNOTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HCells[celli] 			= FGMTable_.interpolateValue4D(HTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2Cells[celli] 			= FGMTable_.interpolateValue4D(H2Tab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2OCells[celli] 		= FGMTable_.interpolateValue4D(H2OTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HO2Cells[celli] 		= FGMTable_.interpolateValue4D(HO2Tab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			H2O2Cells[celli] 		= FGMTable_.interpolateValue4D(H2O2Tab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			OCells[celli] 			= FGMTable_.interpolateValue4D(OTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			O2Cells[celli] 			= FGMTable_.interpolateValue4D(O2Tab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			OHCells[celli] 			= FGMTable_.interpolateValue4D(OHTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			N2Cells[celli] 			= FGMTable_.interpolateValue4D(N2Tab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			NOCells[celli] 			= FGMTable_.interpolateValue4D(NOTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			HRRCells[celli] 		= FGMTable_.interpolateValue4D(HRRTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dCCCells[celli] 		= FGMTable_.interpolateValue4D(dCCTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dCZCells[celli] 		= FGMTable_.interpolateValue4D(dCZTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dZCCells[celli] 		= FGMTable_.interpolateValue4D(dZCTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			dZZCells[celli] 		= FGMTable_.interpolateValue4D(dZZTab,Zi,Ci,ZVi,CVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
		  }

		  
		  if ((Ci == 0) | (Ci ==1) | (Ci >=1.0e0-1e-7) | (Ci <=1.0e-6-1e-7))
		  {
			 sourceYcCells[celli] 	= 0.0;
			 sourceYcVCells[celli] 	= 0.0;
			 sourceH2Cells[celli] 	= 0.0;
			 HRRCells[celli] 		= 0.0;
		  }
	
  }

  // Interpolate for patches
  forAll(Z_.boundaryField(), patchi)
    {
      fvPatchScalarField& pT 			= Temp_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pYcMin 		= YcMin_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pYcMax 		= YcMax_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& psL 			= sL_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& plF 			= lF_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& ptau 			= tau_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pC 			= C_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pdH 			= dH_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pdCC 			= dCC_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pdCZ			= dCZ_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pdZC 			= dZC_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pdZZ 			= dZZ_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pk 		    = k_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pCp 			= Cp_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pW 			= W_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pD 			= D_.boundaryFieldRef()[patchi];
      fvPatchScalarField& pDNO 			= DNO_.boundaryFieldRef()[patchi];
      fvPatchScalarField& psourceYc 	= sourceYc_.boundaryFieldRef()[patchi];
      fvPatchScalarField& psourceYcV 	= sourceYcV_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& psourceH2 	= sourceH2_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& psourceNO 	= sourceNO_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pH 			= H_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pH2 			= H2_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pH2O 			= H2O_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pHO2 			= HO2_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pH2O2 		= H2O2_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pO 			= O_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pO2 			= O2_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pOH 			= OH_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pN2 			= N2_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pNO 			= NO_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pHRR 			= HRR_.boundaryFieldRef()[patchi];
	  fvPatchScalarField& pCV 			= CV_.boundaryFieldRef()[patchi];

      const fvPatchScalarField& pYc 	= Yc_.boundaryField()[patchi];
	  const fvPatchScalarField& pZ 		= Z_.boundaryField()[patchi];
	  const fvPatchScalarField& pYcV 	= YcV_.boundaryField()[patchi];
	  const fvPatchScalarField& pZV 	= ZV_.boundaryField()[patchi];

      forAll(pZ, facei)
      {
		  
		pYcMin[facei] 	= scaleYcTable_.interpolateValue1D(YcMinTab,pZ[facei],ZScalTab);
        pYcMax[facei] 	= scaleYcTable_.interpolateValue1D(YcMaxTab,pZ[facei],ZScalTab);	
        pC[facei] 		= (pYc[facei]-pYcMin[facei])/(pYcMax[facei]-pYcMin[facei]+SMALL);
        pC[facei] 		= max(min(pC[facei], 1.0 ), 0.0 );
		pCV[facei] 		= (pYcV[facei])/(magSqr(pYcMax[facei])+SMALL);
		pCV[facei]	 	= max(min(pCV[facei], 0.25 ), 0.0 );
			
        const double& pZi 	= pZ[facei];
		const double& pCi 	= pC[facei];
		const double& pZVi 	= pZV[facei];
		const double& pCVi 	= pCV[facei];		


		if (pZi < zLean)
		{
			const double ZEff 	= pZi/zLean;
			
			pC[facei]     		= pC[facei]*ZEff;
			pCV[facei]     		= pCV[facei]*ZEff;
			pC[facei] 			= max(min(pC[facei], 1.0 ), 0.0 );
			pCV[facei]	 		= max(min(pCV[facei], 0.25 ), 0.0 );
			
			const double& pCi 	= pC[facei];
			const double& pCVi 	= pCV[facei];
			
			psourceYc[facei] 	= FGMTable_.interpolateValue4D(sourceYcTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceYcV[facei] 	= FGMTable_.interpolateValue4D(sourceYcVTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceH2[facei] 	= FGMTable_.interpolateValue4D(sourceH2Tab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceNO[facei] 	= FGMTable_.interpolateValue4D(sourceNOTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pHRR[facei] 		= FGMTable_.interpolateValue4D(HRRTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pT[facei] 			= FGMTable_.interpolateValue4D(TTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pCp[facei] 			= FGMTable_.interpolateValue4D(CpTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pW[facei] 			= FGMTable_.interpolateValue4D(WTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdH[facei] 			= FGMTable_.interpolateValue4D(dHTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pD[facei] 			= FGMTable_.interpolateValue4D(DTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pDNO[facei] 		= FGMTable_.interpolateValue4D(DNOTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pk[facei] 		    = FGMTable_.interpolateValue4D(kTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH[facei] 			= FGMTable_.interpolateValue4D(HTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2[facei] 			= FGMTable_.interpolateValue4D(H2Tab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2O[facei] 		= FGMTable_.interpolateValue4D(H2OTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2O2[facei] 		= FGMTable_.interpolateValue4D(H2O2Tab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pHO2[facei] 		= FGMTable_.interpolateValue4D(HO2Tab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pO[facei] 			= FGMTable_.interpolateValue4D(OTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pO2[facei] 			= FGMTable_.interpolateValue4D(O2Tab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pOH[facei] 			= FGMTable_.interpolateValue4D(OHTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pN2[facei] 			= FGMTable_.interpolateValue4D(N2Tab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pNO[facei] 			= FGMTable_.interpolateValue4D(NOTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdCC[facei] 		= FGMTable_.interpolateValue4D(dCCTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdCZ[facei] 		= FGMTable_.interpolateValue4D(dCZTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdZC[facei] 		= FGMTable_.interpolateValue4D(dZCTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdZZ[facei] 		= FGMTable_.interpolateValue4D(dZZTab,zLean,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			
			psL[facei] 			= FGMTable_.interpolateValue1D(sLTab,pZi,ZTab);
			plF[facei] 			= FGMTable_.interpolateValue1D(lFTab,pZi,ZTab);
			ptau[facei] 		= FGMTable_.interpolateValue1D(tauTab,pZi,ZTab);
			
			
			pT[facei]     		= FGMTable_.T0_.value() + (pT[facei] - FGMTable_.T0_.value())*ZEff;
			psourceYc[facei]    = 0.0;
			psourceH2[facei]    = psourceH2[facei]*ZEff;
			psourceNO[facei]    = psourceNO[facei]*ZEff;
			psourceYcV[facei]   = 0.0;
			pHRR[facei]    		= pHRR[facei]*ZEff;
			psL[facei]    		= psL[facei]*ZEff;
			plF[facei]    		= plF[facei]*ZEff;
			ptau[facei]    		= ptau[facei]*ZEff;
			pCp[facei]    		= FGMTable_.CpAir_.value() + (pCp[facei] - FGMTable_.CpAir_.value())*ZEff;
			pdH[facei]     		= FGMTable_.dHAir_.value() + (pdH[facei] - FGMTable_.dHAir_.value())*ZEff;
			pW[facei]     		= FGMTable_.WAir_.value() + (pW[facei] - FGMTable_.WAir_.value())*ZEff;
			pD[facei]     		= FGMTable_.DAir_.value() + (pD[facei] - FGMTable_.DAir_.value())*ZEff;
			pDNO[facei]     	= FGMTable_.DNOAir_.value() + (pDNO[facei] - FGMTable_.DNOAir_.value())*ZEff;
			pk[facei]     	    = FGMTable_.kAir_.value() + (pk[facei] - FGMTable_.kAir_.value())*ZEff;
			pH[facei]     		= FGMTable_.HAir_.value() + (pH[facei] - FGMTable_.HAir_.value())*ZEff;
			pH2[facei]     		= FGMTable_.H2Air_.value() + (pH2[facei] - FGMTable_.H2Air_.value())*ZEff;
			pH2O[facei]     	= FGMTable_.H2OAir_.value() + (pH2O[facei] - FGMTable_.H2OAir_.value())*ZEff;
			pH2O2[facei]     	= FGMTable_.H2O2Air_.value() + (pH2O2[facei] - FGMTable_.H2O2Air_.value())*ZEff;
			pHO2[facei]     	= FGMTable_.HO2Air_.value() + (pHO2[facei] - FGMTable_.HO2Air_.value())*ZEff;
			pO[facei]     		= FGMTable_.OAir_.value() + (pO[facei] - FGMTable_.OAir_.value())*ZEff;
			pO2[facei]    	 	= FGMTable_.O2Air_.value() + (pO2[facei] - FGMTable_.O2Air_.value())*ZEff;
			pOH[facei]     		= FGMTable_.OHAir_.value() + (pOH[facei] - FGMTable_.OHAir_.value())*ZEff;
			pN2[facei]     		= FGMTable_.N2Air_.value() + (pN2[facei] - FGMTable_.N2Air_.value())*ZEff;
			pNO[facei]     		= pNO[facei]*ZEff;
			pdCC[facei]     	= FGMTable_.dCCAir_.value() + (pdCC[facei] - FGMTable_.dCCAir_.value())*ZEff;
			pdCZ[facei]     	= FGMTable_.dCZAir_.value() + (pdCZ[facei] - FGMTable_.dCZAir_.value())*ZEff;
			pdZC[facei]     	= FGMTable_.dZCAir_.value() + (pdZC[facei] - FGMTable_.dZCAir_.value())*ZEff;
			pdZZ[facei]     	= FGMTable_.dZZAir_.value() + (pdZZ[facei] - FGMTable_.dZZAir_.value())*ZEff;
		}	
		  
		
		else if (pZi > zRich)
		{
			const double ZEff 	= (pZi-zRich)/(1-zRich);
			
			pC[facei]     		= pC[facei] - pC[facei]*ZEff;
			pCV[facei]     		= pCV[facei] - pCV[facei]*ZEff;
			pC[facei] 			= max(min(pC[facei], 1.0 ), 0.0 );
			pCV[facei] 			= max(min(pCV[facei], 0.25 ), 0.0 );
			
			const double& pCi 	= pC[facei];
			const double& pCVi 	= pCV[facei];
			
			pT[facei] 			= FGMTable_.interpolateValue4D(TTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceNO[facei] 	= FGMTable_.interpolateValue4D(sourceNOTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
		    pCp[facei] 			= FGMTable_.interpolateValue4D(CpTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pW[facei] 			= FGMTable_.interpolateValue4D(WTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdH[facei] 			= FGMTable_.interpolateValue4D(dHTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pD[facei] 			= FGMTable_.interpolateValue4D(DTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pDNO[facei] 	    = FGMTable_.interpolateValue4D(DNOTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pk[facei] 		    = FGMTable_.interpolateValue4D(kTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH[facei] 			= FGMTable_.interpolateValue4D(HTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2[facei] 			= FGMTable_.interpolateValue4D(H2Tab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2O[facei] 		= FGMTable_.interpolateValue4D(H2OTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2O2[facei] 		= FGMTable_.interpolateValue4D(H2O2Tab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pHO2[facei] 		= FGMTable_.interpolateValue4D(HO2Tab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pO[facei] 			= FGMTable_.interpolateValue4D(OTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pO2[facei] 			= FGMTable_.interpolateValue4D(O2Tab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pOH[facei] 			= FGMTable_.interpolateValue4D(OHTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pN2[facei] 			= FGMTable_.interpolateValue4D(N2Tab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pNO[facei] 			= FGMTable_.interpolateValue4D(NOTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdCC[facei] 		= FGMTable_.interpolateValue4D(dCCTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdCZ[facei] 		= FGMTable_.interpolateValue4D(dCZTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdZC[facei] 		= FGMTable_.interpolateValue4D(dZCTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdZZ[facei] 		= FGMTable_.interpolateValue4D(dZZTab,zRich,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			
			psourceYc[facei]    = 0.0;
			psourceYcV[facei]   = 0.0;
			psourceH2[facei]    = 0.0;
			pHRR[facei]         = 0.0;
			psL[facei]     		= 0.0;
			plF[facei]     		= 0.0;
			ptau[facei]     	= 0.0;
			
			pT[facei]    		= pT[facei]  + (FGMTable_.T0_.value() - pT[facei])*ZEff;
			psourceNO[facei]    = psourceNO[facei]  - psourceNO[facei]*ZEff;
			pCp[facei]     		= pCp[facei]  + (FGMTable_.CpFuel_.value() - pCp[facei])*ZEff;
			pdH[facei]     		= pdH[facei]  + (FGMTable_.dHFuel_.value() - pdH[facei])*ZEff;
			pW[facei]     		= pW[facei]  + (FGMTable_.WFuel_.value() - pW[facei])*ZEff;
			pD[facei]     		= pD[facei]  + (FGMTable_.DFuel_.value() - pD[facei])*ZEff;
			pDNO[facei]     	= pDNO[facei]  + (FGMTable_.DNOFuel_.value() - pD[facei])*ZEff;
			pk[facei]     	    = pk[facei]  + (FGMTable_.kFuel_.value() - pk[facei])*ZEff;
			pH[facei]    	    = pH[facei]  + (FGMTable_.HFuel_.value() - pH[facei])*ZEff;
			pH2[facei]     		= pH2[facei]  + (FGMTable_.H2Fuel_.value() - pH2[facei])*ZEff;
			pH2O[facei]     	= pH2O[facei]  + (FGMTable_.H2OFuel_.value() - pH2O[facei])*ZEff;
			pH2O2[facei]     	= pH2O2[facei]  + (FGMTable_.H2O2Fuel_.value() - pH2O2[facei])*ZEff;
			pHO2[facei]     	= pHO2[facei]  + (FGMTable_.HO2Fuel_.value() - pHO2[facei])*ZEff;
			pO[facei]     		= pO[facei]  + (FGMTable_.OFuel_.value() - pO[facei])*ZEff;
			pO2[facei]     		= pO2[facei]  + (FGMTable_.O2Fuel_.value() - pO2[facei])*ZEff;
			pOH[facei]     		= pOH[facei]  + (FGMTable_.OHFuel_.value() - pOH[facei])*ZEff;
			pN2[facei]     		= pN2[facei]  + (FGMTable_.N2Fuel_.value() - pN2[facei])*ZEff;
			pNO[facei]     		= pNO[facei]  + - pNO[facei]*ZEff;
			pdCC[facei]     	= pdCC[facei]  + (FGMTable_.dCCFuel_.value() - pdCC[facei])*ZEff;
			pdCZ[facei]     	= pdCZ[facei]  + (FGMTable_.dCZFuel_.value() - pdCZ[facei])*ZEff;
			pdZC[facei]     	= pdZC[facei]  + (FGMTable_.dZCFuel_.value() - pdZC[facei])*ZEff;
			pdZZ[facei]     	= pdZZ[facei]  + (FGMTable_.dZZFuel_.value() - pdZZ[facei])*ZEff;	  					  
		}

        else
		{
			psourceYc[facei] 	= FGMTable_.interpolateValue4D(sourceYcTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceYcV[facei] 	= FGMTable_.interpolateValue4D(sourceYcVTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceH2[facei] 	= FGMTable_.interpolateValue4D(sourceH2Tab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceNO[facei] 	= FGMTable_.interpolateValue4D(sourceNOTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			psourceYc[facei] 	= max(psourceYc[facei], 0.0);
			psourceYcV[facei] 	= max(psourceYcV[facei], 0.0);

			psL[facei] 			= FGMTable_.interpolateValue1D(sLTab,pZi,ZTab);
			plF[facei] 			= FGMTable_.interpolateValue1D(lFTab,pZi,ZTab);
			ptau[facei] 		= FGMTable_.interpolateValue1D(tauTab,pZi,ZTab);
					  
			pT[facei] 			= FGMTable_.interpolateValue4D(TTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pk[facei] 		    = FGMTable_.interpolateValue4D(kTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pCp[facei] 			= FGMTable_.interpolateValue4D(CpTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdH[facei] 			= FGMTable_.interpolateValue4D(dHTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pW[facei] 			= FGMTable_.interpolateValue4D(WTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pD[facei] 			= FGMTable_.interpolateValue4D(DTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pDNO[facei] 		= FGMTable_.interpolateValue4D(DNOTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH[facei] 			= FGMTable_.interpolateValue4D(HTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2[facei] 			= FGMTable_.interpolateValue4D(H2Tab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2O[facei] 		= FGMTable_.interpolateValue4D(H2OTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pH2O2[facei] 		= FGMTable_.interpolateValue4D(H2O2Tab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pHO2[facei] 		= FGMTable_.interpolateValue4D(HO2Tab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pO[facei] 			= FGMTable_.interpolateValue4D(OTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pO2[facei] 			= FGMTable_.interpolateValue4D(O2Tab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pOH[facei] 			= FGMTable_.interpolateValue4D(OHTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pN2[facei] 			= FGMTable_.interpolateValue4D(N2Tab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pNO[facei] 			= FGMTable_.interpolateValue4D(NOTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pHRR[facei] 		= FGMTable_.interpolateValue4D(HRRTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdCC[facei] 		= FGMTable_.interpolateValue4D(dCCTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdCZ[facei] 		= FGMTable_.interpolateValue4D(dCZTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdZC[facei] 		= FGMTable_.interpolateValue4D(dZCTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
			pdZZ[facei] 		= FGMTable_.interpolateValue4D(dZZTab,pZi,pCi,pZVi,pCVi,ZTab,CTab,ZVTab,CVTab,nPointsC,nPointsZv,nPointsCv);
															
		}

       	if ((pCi == 0) | (pCi ==1) | (pCi >=1.0e0-1e-7) | (pCi <=1.0e-6-1e-7))
	    {
		  psourceYc[facei] 		= 0.0;
		  psourceYcV[facei] 	= 0.0;
		  psourceH2[facei] 		= 0.0;
		  pHRR[facei] 			= 0.0;
	    }
      }
    }
	
//updating Cv
Cv_ == Cp_ - R0/W_;

//this->T_ = T0 + (he_ - dH_)/Cp_;
this->T_ = Temp_;
this->psi_ = W_/(R0*Temp_);

Info << "min/max T: "<<min(T_).value()<<"/"<<max(T_).value() << endl;

this->mu_ = this->mu(); //As*sqrt(this->T_)/(1.0+Ts/this->T_)
this->alpha_ = this->alphahe(); // this->Cv()*this->mu_*(1.32+1.77*RR/(Wmix_*this->Cv()))/Cp_;
}

tmp<volScalarField> FGMThermo::delta() const
{
    const fvMesh& mesh = this->T_.mesh();
    const scalarField& VCells = mesh.V();

    tmp<volScalarField> tDelta
    (
        new volScalarField
        (
            IOobject
            (
                "delta",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("delta",dimLength,1.0)
        )
    );

    volScalarField& delta = tDelta.ref();
    scalarField& deltaCells = delta.ref();

    forAll(deltaCells, celli)
        {
            deltaCells[celli] = pow(VCells[celli],1.0/3.0);
        }

   forAll(delta.boundaryField(), patchi)
        {
            scalarField& deltap = delta.boundaryFieldRef()[patchi];

            forAll(deltap, facei)
                {
                    label celli = mesh.boundary()[patchi].faceCells()[facei];
                    deltap[facei] = delta[celli]; //1.0;
                }
	}

    return tDelta;
}

tmp<volScalarField> FGMThermo::Da(const volScalarField& k, const volScalarField& delta) const
{
	
    const fvMesh& mesh = this->T_.mesh();
    tmp<volScalarField> tDa
    (
        new volScalarField
        (
            IOobject
            (
                "Da",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );
    volScalarField& Da = tDa.ref();
    const dimensionedScalar smallDim("smallDim",dimensionSet(0,2,-1,0,0,0,0),1.0E-08);
    Da = (delta*sL_)/(lF_*sqrt(2.0*k/3.0)+smallDim);
    return tDa;
}

tmp<volScalarField> FGMThermo::Ka(const volScalarField& k, const volScalarField& delta) const
{

    const fvMesh& mesh = this->T_.mesh();
    tmp<volScalarField> tKa
    (
        new volScalarField
        (
            IOobject
            (
                "Ka",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );
    volScalarField& Ka = tKa.ref();
    Ka = pow(sqrt(2.0*k/3.0)/sL_,1.5)/(sqrt(delta/lF_)+1.0E-08);
    return tKa;
}

tmp<volScalarField> FGMThermo::C2(const volScalarField& k, const volScalarField& T, const volScalarField& delta)
{
   const scalar theta5 = 0.75;
   const scalarList& ZTab = FGMTable_.Z_table;
   const scalarList& sLTab = FGMTable_.sL_table;
   const scalarList& lFTab = FGMTable_.lF_table;
   const scalarList& tauTab = FGMTable_.tau_table;
   const dimensionedScalar T0(FGMTable_.lookup("T0"));
   const dimensionedScalar zLean(FGMTable_.lookup("zLean"));

   const fvMesh& mesh = this->T_.mesh();

   tmp<volScalarField> tC2
    (
        new volScalarField
        (
            IOobject
            (
                "C2",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0,0,-1,0,0,0,0)
        )
    );

    volScalarField& C2 = tC2.ref(); 
	  
	volScalarField tau
   (
        IOobject
        (
            "tau",
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0)
    ); 
	
	
    const scalarField& ZCells = Z_.internalField();
    const scalarField& kCells = k.internalField();
	//const scalarField& TCells = T_.internalField();
	scalarField& sLi = sL_.primitiveFieldRef();
	scalarField& lFi = lF_.primitiveFieldRef();
	scalarField& taui = tau_.primitiveFieldRef();
    const scalarField& delCells = delta.internalField();

    scalarField& C2Cells = C2.ref();
    
    forAll(ZCells,celli)
    {
	  const scalar& Zi =  ZCells[celli];
	  const scalar& ki =  kCells[celli];
	  const scalar& deli = delCells[celli];
	  
	  if (Zi<zLean.value())
          {
             C2Cells[celli] = sqrt(2.0*ki/3.0)/deli;
          }
		  
	  else
          {
			 sLi[celli] = FGMTable_.interpolateValue1D(sLTab,Zi,ZTab);
			 lFi[celli] = FGMTable_.interpolateValue1D(lFTab,Zi,ZTab);
			 taui[celli] = FGMTable_.interpolateValue1D(tauTab,Zi,ZTab);
             scalar Kai = pow(sqrt(2.0*ki/3.0)/sLi[celli],1.5)/(sqrt(deli/lFi[celli])+1.0E-08);
             scalar Dai = (deli*sLi[celli])/(lFi[celli]*sqrt(2.0*ki/3.0)+1.0E-08);

             C2Cells[celli] 
			 = 
			   (1.0 - exp(-theta5*deli/lFi[celli]))*(2.0*0.79*taui[celli]*sLi[celli]/lFi[celli]
             + ( 1.5*sqrt(Kai)/(1.0+sqrt(Kai))     //C3
             -  taui[celli]*Dai*1.1/(pow(1.0+Kai,0.4)))*(2.0/3.0*sqrt(2.0*ki/3.0)/deli));  
          }
       }
	   
	forAll(Z_.boundaryField(), patchi)
       {
            const fvPatchScalarField& ZBound = Z_.boundaryField()[patchi];
            fvPatchScalarField& C2Bound = C2.boundaryFieldRef()[patchi];

            forAll(ZBound, facei)
            {
                label celli = mesh.boundary()[patchi].faceCells()[facei];
                C2Bound[facei] = C2Cells[celli];
            }

       }
	   return tC2;
}
	 
tmp<volScalarField> FGMThermo::zvSGSdiss(const volScalarField& muT, const volScalarField& delta)
{

   const fvMesh& mesh = this->T_.mesh();

   tmp<volScalarField> tZvDiss
    (
        new volScalarField
        (
            IOobject
            (
                "ZvDiss",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(1,-3,-1,0,0,0,0)
        )
    );

    volScalarField& ZvDiss = tZvDiss.ref();

    const scalarField& muCells = muT.internalField();
    const scalarField& delCells = delta.internalField();


    scalarField& ZvDissCells = ZvDiss.ref();

    forAll(ZvDissCells,celli)
    {
       const scalar& mui =  muCells[celli];
       const scalar& deli = delCells[celli];

       ZvDissCells[celli] = mui/(deli*deli);
    }

    forAll(T_.boundaryField(), patchi)
    {

         const fvPatchScalarField& TBound = T_.boundaryField()[patchi];
         fvPatchScalarField& ZvDissBound = ZvDiss.boundaryFieldRef()[patchi];

         forAll(TBound, facei)
         {
             label celli = mesh.boundary()[patchi].faceCells()[facei];
             ZvDissBound[facei] = ZvDissCells[celli];
         }

    }

   return tZvDiss;
}	

tmp<volScalarField> FGMThermo::filter(const tmp<volScalarField>& phi)
{

     const fvMesh& mesh = this->T_.mesh();
     tmp<volScalarField> tphiFilt
     (
        new volScalarField
        (
           IOobject
           (
               "phiFilt",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           phi,
           Yc_.boundaryField().types()
        )

     );


     volScalarField& phiFilt = tphiFilt.ref();

     tmp<volScalarField> tfilterCoeff
     (
       new volScalarField
       (
        IOobject
        (
           "filterCoeff",
           mesh.time().timeName(),
           mesh
        ),
        mesh,
        dimensionedScalar("filterCoeff", dimLength*dimLength/dimTime, 0.0),
        calculatedFvPatchScalarField::typeName
       )
     );

     volScalarField& filterCoeff = tfilterCoeff.ref();
     dimensionedScalar dt = mesh.time().deltaT();
     const dimensionedScalar dt1("dt1",dimTime,1.0);

     filterCoeff.ref() = pow(mesh.V(), 2.0/3.0)/(6.0*dt);

     const scalarField& phiCells = phiFilt.internalField();
     const scalarField& coeffCells = filterCoeff.internalField();

     forAll(Yc_.boundaryField(), patchi)
     {

          const fvPatchScalarField& YcBound = Yc_.boundaryField()[patchi];
          fvPatchScalarField& phiBound = phiFilt.boundaryFieldRef()[patchi];
          fvPatchScalarField& coeffBound = filterCoeff.boundaryFieldRef()[patchi];


          forAll(YcBound, facei)
          {
              label celli = mesh.boundary()[patchi].faceCells()[facei];
              phiBound[facei] = phiCells[celli];
              coeffBound[facei] = coeffCells[celli];
          }

     }


     solve
     (
          fvm::ddt(phiFilt)
        - filterCoeff*fvm::laplacian(phiFilt)
     );

     phi.clear();

     return tphiFilt;
}

tmp<volScalarField> FGMThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> the
    (
        new volScalarField
        (
            IOobject
            (
                "he",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& he = the.ref();
    scalarField& heCells = he.ref();
    //const scalarField& pCells = p.internalField();
    const scalarField& TCells = T.internalField();
    const dimensionedScalar T0(FGMTable_.lookup("T0"));

    forAll(heCells, celli)
    {
        heCells[celli] = dH_.internalField()[celli] +
                         Cp_.internalField()[celli]*(TCells[celli] - T0.value());

    }

    forAll(he.boundaryField(), patchi)
    {
        scalarField& hep = he.boundaryFieldRef()[patchi];
        //const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(hep, facei)
        {
            hep[facei] = dH_.boundaryField()[patchi][facei]
                 + Cp_.boundaryField()[patchi][facei]*(Tp[facei]-T0.value());
        }
    }

    return the;
}

tmp<scalarField> FGMThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th.ref();
    const dimensionedScalar T0(FGMTable_.lookup("T0"));
    forAll(T, celli)
    {   
        h[celli] = dH_[cells[celli]] + Cp_[cells[celli]]*(T[celli] - T0.value());   
    }

    return th;
}

tmp<scalarField> FGMThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th.ref();
    const dimensionedScalar T0(FGMTable_.lookup("T0"));

    forAll(T, facei)
    {
        h[facei] = dH_.boundaryField()[patchi][facei]
                 + Cp_.boundaryField()[patchi][facei]*(T[facei]-T0.value());
    }

    return th;
}

tmp<volScalarField> FGMThermo::hs() const
{
    const dimensionedScalar T0(FGMTable_.lookup("T0"));
    return (this->Cp_*(this->T_ - T0));
}

tmp<volScalarField> FGMThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
	    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> ths
    (
        new volScalarField
        (
            IOobject
            (
                "hs",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& hs = ths.ref();
    scalarField& hsCells = hs.ref();
    //const scalarField& pCells = p.internalField();
    const scalarField& TCells = T.internalField();
    const dimensionedScalar T0(FGMTable_.lookup("T0"));

    forAll(hsCells, celli)
    {
        hsCells[celli] = Cp_.internalField()[celli]*(TCells[celli] - T0.value());

    }

    forAll(hs.boundaryField(), patchi)
    {
        scalarField& hsp = hs.boundaryFieldRef()[patchi];
        //const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tp = T.boundaryField()[patchi];

        forAll(hsp, facei)
        {
            hsp[facei] = Cp_.boundaryField()[patchi][facei]*(Tp[facei]-T0.value());
        }
    }

    return ths;
}

tmp<scalarField> FGMThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th.ref();
    const dimensionedScalar T0(FGMTable_.lookup("T0"));
    forAll(T, celli)
    {   
        h[celli] = Cp_[cells[celli]]*(T[celli] - T0.value());   
    }

    return th;
}


tmp<scalarField> FGMThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th.ref();
    const dimensionedScalar T0(FGMTable_.lookup("T0"));


    forAll(T, facei)
    {   
        h[facei] = Cp_.boundaryField()[patchi][facei]*(T[facei]-T0.value());
    }

    return th;    
}

tmp<Foam::volScalarField> FGMThermo::ha() const
{
	return this->he_;
}

tmp<volScalarField> FGMThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
	return he(p,T);
}


tmp<scalarField> FGMThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return he(T,cells);
}

tmp<scalarField> FGMThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return he(T,patchi);
}


tmp<volScalarField> FGMThermo::hc() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& hcf = thc.ref();
    scalarField& hcCells = hcf.ref();

    forAll(hcCells, celli)
    {
        hcCells[celli] = dH_.internalField()[celli];
    }

    forAll(hcf.boundaryField(), patchi)
    {
        scalarField& hcp = hcf.boundaryFieldRef()[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = dH_.boundaryField()[patchi][facei];
        }
    }

    return thc;
}


tmp<scalarField> FGMThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& Cp = tCp.ref();

    forAll(T, facei)
    {   
        Cp[facei] = Cp_.boundaryField()[patchi][facei];
    }

    return tCp;
}


tmp<scalarField> FGMThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& Cv = tCv.ref();
    const dimensionedScalar R0(FGMTable_.lookup("R0"));

    forAll(T, facei)
    {
        Cv[facei] = Cp_.boundaryField()[patchi][facei]
                  - R0.value()/W_.boundaryField()[patchi][facei];
    }

    return tCv;
}


tmp<scalarField> FGMThermo::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tgamma(new scalarField(T.size()));
    scalarField& gamma = tgamma.ref();

    gamma = Cp(T,patchi)/Cv(T,patchi);

    return tgamma;
}


tmp<volScalarField> FGMThermo::gamma() const
{
    return volScalarField::New("gamma", Cp_/Cv_);
}


tmp<scalarField> FGMThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpv(new scalarField(T.size()));
    scalarField& cpv = tCpv.ref();

    forAll(T, facei)
    {
        cpv[facei] =
            this->Cp_.boundaryField()[patchi][facei];
    }

    return tCpv;
}


tmp<volScalarField> FGMThermo::Cpv() const
{
        return Cp_;
}


tmp<volScalarField> FGMThermo::THE
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
	return (T0 + (h - dH_)/Cp_);
}


tmp<scalarField> FGMThermo::THE
(
    const scalarField& h,
    const scalarField& T0,
    const labelList& cells
) const
{
    tmp<scalarField> tT(new scalarField(h.size()));
    scalarField& T = tT.ref();

    forAll(h, celli)
    {
        T[celli] = T0[celli] + (h[celli]-dH_[cells[celli]])/Cp_[cells[celli]];
    }

    return tT;
}


tmp<scalarField> FGMThermo::THE
(
    const scalarField& h,
    const scalarField& T0,
    const label patchi
) const
{
    tmp<scalarField> tT(new scalarField(h.size()));
    scalarField& T = tT.ref();
    forAll(h, facei)
    {
        T[facei] = T0[facei] + (h[facei]-dH_.boundaryField()[patchi][facei])
                               /Cp_.boundaryField()[patchi][facei];
    }

    return tT;
}


tmp<volScalarField> FGMThermo::W() const
{
    return this->W_;
}


tmp<scalarField> FGMThermo::Temp
(
    const label patchi
) const
{
    return this->Temp_.boundaryField()[patchi];
}

tmp<volScalarField> FGMThermo::Temp() const
{
    return this->Temp_;
}


tmp<scalarField> FGMThermo::W
(
    const label patchi
) const
{
    return this->W_.boundaryField()[patchi];
}



tmp<volScalarField> FGMThermo::kappa() const
{
    return volScalarField::New("kappa", Cp_*this->alpha_);
}

tmp<scalarField> FGMThermo::kappa
(
    const label patchi
) const
{
    return
        Cp
        (
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}

tmp<volScalarField> FGMThermo::mu() const
{
    const fvMesh& mesh = this->T_.mesh();
	const dimensionedScalar As("As",dimensionSet(1,-1,-1,-0.5,0,0,0),1.67212e-6);
    const dimensionedScalar Ts("Ts",dimensionSet(0,0,0,1,0,0,0),170.672);

    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(1,-1,-1,0,0,0,0)
        )
    );
    volScalarField& mu = tmu.ref();
	mu = As*sqrt(this->T_)/(1.0+Ts/this->T_);

    return tmu;
}

tmp<volScalarField> FGMThermo::alphahe() const
{
    const fvMesh& mesh = this->T_.mesh();
    const dimensionedScalar R0(FGMTable_.lookup("R0"));

    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                "alpha",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(1,-1,-1,0,0,0,0)
        )
    );

    volScalarField& alpha = talpha.ref();
    alpha = Cv_*mu()*(1.32+1.77*R0/(W_*Cv_))/Cp_;
    return talpha;
}


tmp<scalarField> FGMThermo::alphahe(const label patchi) const
{
    return this->alpha_.boundaryField()[patchi];
}


tmp<volScalarField> FGMThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New("kappaEff", Cp_*(this->alpha_ + alphat));
}


tmp<scalarField> FGMThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        Cp
        (
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(this->alpha_.boundaryField()[patchi] + alphat);
}

tmp<volScalarField> FGMThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New("alphaEff", this->alpha_ + alphat);
}


tmp<scalarField> FGMThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return (this->alpha_.boundaryField()[patchi] + alphat);
}


bool FGMThermo::incompressible() const
        {
            //Switch icprs(combustionProperties.lookup("LowMach"));
            return true;
        }

bool FGMThermo::read()
{
    if (basicThermo::implementation::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


} //End namespace Foam
// ************************************************************************* //
