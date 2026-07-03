/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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
       
#include "scaleYc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

scaleYc::scaleYc
(
 const fvMesh& mesh
)
  :
  IOdictionary
  (
    IOobject
    (
     "scalingYcTable",
     mesh.time().constant(),
     mesh,
     IOobject::MUST_READ_IF_MODIFIED,
     IOobject::NO_WRITE
     )
  ),
  mesh_(mesh),
  
  //Look-up Variables in the table
  ZScal(lookup("ZScal")),
  YcMax(lookup("YcMax")),
  YcMin(lookup("YcMin"))

{
  Info << "\nReading scalingYcTable" << nl << endl;
  Info << "Table length: " << ZScal.size() << endl;
  Info << "Table contents:" << endl;
  Info << "{" << endl;
  Info << "Mixture Fraction: Z" << endl;
  Info << "Max Yc(Z): YcMax" << endl;
  Info << "Min Yc(Z): YcMin" << endl;
  Info << "}" << endl;
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

scaleYc::~scaleYc()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::scaleYc::interpolateValue1D
(
    const List<scalar>& table,   // length = Nz
    const scalar zValue,         // query value
    const List<scalar>& Zvec     // coordinate vector
) const
{
    scalar smallValue = 1e-12;

    const label Nz = table.size();

    if (Nz != Zvec.size())
    {
        FatalErrorInFunction
            << "table and Zvec must have the same length" << nl
            << "table.size() = " << table.size()
            << ", Zvec.size() = " << Zvec.size() << nl
            << exit(FatalError);
    }

    if (Nz < 2)
    {
        FatalErrorInFunction
            << "Need at least two points for interpolation" << nl
            << exit(FatalError);
    }

    // Clamp query value inside the table range
    scalar zq = zValue;
    if (zq <= Zvec[0]) return table[0];
    if (zq >= Zvec[Nz-1]) return table[Nz-1];

    // Find lower index i such that Zvec[i] <= zq <= Zvec[i+1]
    label i = 0;
    for (label k = 1; k < Nz; ++k)
    {
        if (zq <= Zvec[k])
        {
            i = k - 1;
            break;
        }
    }

    // Linear interpolation
    scalar z0 = Zvec[i];
    scalar z1 = Zvec[i+1];
    scalar f0 = table[i];
    scalar f1 = table[i+1];

    scalar t = (zq - z0) / max(z1 - z0, smallValue);

    return f0 + t * (f1 - f0);
}