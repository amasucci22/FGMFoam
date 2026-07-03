/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
 
Class
    Foam::lookupFGM
 
Description
    Class to look up variables from FGMTable and perform the 4D 
	interpolation.

    @author:  Antonio Masucci
    @created: 01/05/2026
 
SourceFiles
    lookupFGM.C
 
\*---------------------------------------------------------------------------*/

#ifndef lookupFGM_H
#define lookupFGM_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalar.H"
#include "scalarField.H"
#include "volFields.H"
#include "PtrList.H"
#include "zeroGradientFvPatchFields.H"
#include "OFstream.H"
#include "IFstream.H"

#include "List.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class lookupFGM
  :
    public IOdictionary
{
  
protected:

    // Protected data
  
        //- Reference to the mesh database
        const fvMesh& mesh_; 
private:
    
  
public:
     //- Return the instantiated type name
      static word typeName()
       {
         return "lookupFGM";
       }
  
    // Constructors

    //- Construct from components
    lookupFGM(const fvMesh&);

  //- Destructor
    ~lookupFGM();

    // Access general FGM table data
	//Air and fuel values
	const dimensionedScalar zLean_;
	const dimensionedScalar zRich_;
	const dimensionedScalar T0_;
	const dimensionedScalar R0_;
	const dimensionedScalar pRef_;
	
	const dimensionedScalar CpAir_;
	const dimensionedScalar DAir_;
	const dimensionedScalar DNOAir_;
	const dimensionedScalar WAir_;
	const dimensionedScalar kAir_;
	const dimensionedScalar dHAir_;
	const dimensionedScalar HAir_;
	const dimensionedScalar H2Air_;
	const dimensionedScalar H2OAir_;
	const dimensionedScalar H2O2Air_;
	const dimensionedScalar HO2Air_;
	const dimensionedScalar OAir_;
	const dimensionedScalar O2Air_;
	const dimensionedScalar OHAir_;
	const dimensionedScalar N2Air_;
	const dimensionedScalar OHEXAir_;
	const dimensionedScalar dCCAir_;
	const dimensionedScalar dCZAir_;
	const dimensionedScalar dZCAir_;
	const dimensionedScalar dZZAir_;
	
	const dimensionedScalar CpFuel_;
	const dimensionedScalar DFuel_;
	const dimensionedScalar DNOFuel_;
	const dimensionedScalar WFuel_;
	const dimensionedScalar kFuel_;
	const dimensionedScalar dHFuel_;
	const dimensionedScalar HFuel_;
	const dimensionedScalar H2Fuel_;
	const dimensionedScalar H2OFuel_;
	const dimensionedScalar H2O2Fuel_;
	const dimensionedScalar HO2Fuel_;
	const dimensionedScalar OFuel_;
	const dimensionedScalar O2Fuel_;
	const dimensionedScalar OHFuel_;
	const dimensionedScalar N2Fuel_;
	const dimensionedScalar OHEXFuel_;
	const dimensionedScalar dCCFuel_;
	const dimensionedScalar dCZFuel_;
	const dimensionedScalar dZCFuel_;
	const dimensionedScalar dZZFuel_;
	
	//Table values
	const scalarList Z_table;
    const scalarList C_table;
	const scalarList ZV_table;
    const scalarList CV_table;
	const scalarList sL_table;
	const scalarList lF_table;
	const scalarList tau_table;
    const scalarList sourceYc_table;
	const scalarList sourceYcV_table;
	const scalarList sourceH2_table;
	const scalarList sourceNO_table;
	const scalarList dCC_table;
    const scalarList dCZ_table;
	const scalarList dZC_table;
	const scalarList dZZ_table;
	const scalarList k_table;
	const scalarList dH_table;
    const scalarList Cp_table;
	const scalarList W_table;
	const scalarList rho_table;
	const scalarList Temp_table;
	const scalarList D_table;  
	const scalarList DNO_table;  
	const scalarList H_table;
	const scalarList H2_table;
	const scalarList H2O_table;
	const scalarList H2O2_table;
	const scalarList HO2_table;
	const scalarList O_table;
	const scalarList O2_table;
	const scalarList OH_table;
	const scalarList N2_table;
	const scalarList NO_table;
	const scalarList OHEX_table;
	const scalarList HRR_table;
	label nZ_;
	label nC_;
	label nZv_;
    label nCv_;
	
  
    //- Return const access to the mesh database
    inline const fvMesh& mesh() const;
  
    //- Interpolate value according FGM table
      scalar interpolateValue4D
        (
	    const List<scalar>&,
            scalar,
			scalar,
			scalar,
			scalar,
			const List<scalar>&,
            const List<scalar>&,
			const List<scalar>&,
            const List<scalar>&,
			const int& ,
			const int& ,
			const int& 
        ) const;  
		
		scalar interpolateValue1D
        (
	    const List<scalar>&,
			scalar,
			const List<scalar>&
        ) const;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
