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
    Foam::scaleYc
 
Description
    Class to look up Yc.b to create the scaled progress variable C \in [0,1]

    @author:  Antonio Masucci
    @created: 01/05/2026
 
SourceFiles
    scaleYc.C
 
\*---------------------------------------------------------------------------*/

#ifndef scaleYc_H
#define scaleYc_H

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
class scaleYc
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
         return "scaleYc";
       }
  
    // Constructors

    //- Construct from components
    scaleYc(const fvMesh&);

  //- Destructor
    ~scaleYc();

    // Access general FGM table data
    List<scalar> ZScal;
	List<scalar> YcMax;
	List<scalar> YcMin;

	
  
    //- Return const access to the mesh database
    inline const fvMesh& mesh() const;
  
    //- Interpolate value according FGM table
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
