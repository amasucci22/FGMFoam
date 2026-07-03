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
       
#include "lookupFGM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lookupFGM::lookupFGM
(
 const fvMesh& mesh
)
  :
  IOdictionary
  (
    IOobject
    (
     "fgmProperties",
     mesh.time().constant(),
     mesh,
     IOobject::MUST_READ_IF_MODIFIED,
     IOobject::NO_WRITE
     )
  ),
  mesh_(mesh),
  
  //Look-up Variables in the table
  //Default Values in Air and Fuel
  zLean_(lookup("zLean")),
  zRich_(lookup("zRich")),
  T0_(lookup("T0")),
  R0_(lookup("R0")),
  pRef_(lookup("pRef")),
  
  CpAir_(lookup("CpAir")),
  DAir_(lookup("DAir")),
  DNOAir_(lookup("DNOAir")),						
  WAir_(lookup("WAir")),
  kAir_(lookup("kAir")),
  dHAir_(lookup("dHAir")),
  HAir_(lookup("YHAir")),
  H2Air_(lookup("YH2Air")),
  H2OAir_(lookup("YH2OAir")),
  H2O2Air_(lookup("YH2O2Air")),
  HO2Air_(lookup("YHO2Air")),
  OAir_(lookup("YOAir")),
  O2Air_(lookup("YO2Air")),
  OHAir_(lookup("YOHAir")),
  N2Air_(lookup("YN2Air")),
  dCCAir_(lookup("dCCAir")),
  dCZAir_(lookup("dCZAir")),
  dZCAir_(lookup("dZCAir")),
  dZZAir_(lookup("dZZAir")),
  
  CpFuel_(lookup("CpFuel")),
  DFuel_(lookup("DFuel")),
  DNOFuel_(lookup("DNOFuel")),						  
  WFuel_(lookup("WFuel")),
  kFuel_(lookup("kFuel")),
  dHFuel_(lookup("dHFuel")),
  HFuel_(lookup("YHFuel")),
  H2Fuel_(lookup("YH2Fuel")),
  H2OFuel_(lookup("YH2OFuel")),
  H2O2Fuel_(lookup("YH2O2Fuel")),
  HO2Fuel_(lookup("YHO2Fuel")),
  OFuel_(lookup("YOFuel")),
  O2Fuel_(lookup("YO2Fuel")),
  OHFuel_(lookup("YOHFuel")),
  N2Fuel_(lookup("YN2Fuel")),
  dCCFuel_(lookup("dCCFuel")),
  dCZFuel_(lookup("dCZFuel")),
  dZCFuel_(lookup("dZCFuel")),
  dZZFuel_(lookup("dZZFuel")),
  
  //Variables
  Z_table(lookup("Z")),
  C_table(lookup("C")),
  ZV_table(lookup("Zv")),
  CV_table(lookup("Cv")),
  sL_table(lookup("sL")),
  lF_table(lookup("lF")),
  tau_table(lookup("tau")),
  sourceYc_table(lookup("omegaYc")),
  sourceYcV_table(lookup("omegaYcV")),
  sourceH2_table(lookup("omegaH2")),
  sourceNO_table(lookup("omegaNO")),							
  dCC_table(lookup("dYcYc")),
  dCZ_table(lookup("dYcZ")),
  dZC_table(lookup("dZYc") ),
  dZZ_table(lookup("dZZ")),
  k_table(lookup("k")),
  dH_table(lookup("dH")),
  Cp_table(lookup("Cp")),
  W_table(lookup("W")),
  rho_table(lookup("rho")),
  Temp_table(lookup("T")),
  D_table(lookup("D")),
  DNO_table(lookup("DNO")),					   
  H_table(lookup("Yk_H")),
  H2_table(lookup("Yk_H2")),
  H2O_table(lookup("Yk_H2O")),
  H2O2_table(lookup("Yk_H2O2")),
  HO2_table(lookup("Yk_HO2")),
  O_table(lookup("Yk_O")),
  O2_table(lookup("Yk_O2")),
  OH_table(lookup("Yk_OH")),
  N2_table(lookup("Yk_N2")),
  NO_table(lookup("Yk_NO")),							
  HRR_table(lookup("HRR")),
  nZ_(readLabel(lookup("nZ"))),
  nC_(readLabel(lookup("nC"))),
  nZv_(readLabel(lookup("nZv"))),
  nCv_(readLabel(lookup("nCv")))

{
  Info << "\nLook-up Table initialization" << endl;
  Info << "Reading fgmProperties" << endl;
  Info << "Table length: " << sourceYc_table.size() << endl;
  Info << "Table contents:" << endl;
  Info << "{" << endl;
  Info << "Progress Variable: C" << endl;
  Info << "Mixture Fraction: Z" << endl;
  Info << "Progress Variable Var: CV" << endl;
  Info << "Mixture Fraction Var: ZV" << endl;
  Info << "Laminar Flame Speed: sL" << endl;
  Info << "Laminar Flame Thickness: lF" << endl;
  Info << "Heat Release Factor: tau" << endl;
  Info << "Source Term of Yc: sourceYc" << endl;
  Info << "Source Term of YcV: sourceYcV" << endl;
  Info << "Source Term of H2: sourceH2" << endl;
  Info << "Species Mass fractions Yk: Yk" << endl;
  Info << "Thermal Conductivity: k" << endl;
  Info << "Mean Molecular Weight: W" << endl;
  Info << "Temperature: T" << endl;
  Info << "Density: rho" << endl;
  Info << "Heat Release Rate: HRR" << endl;
  Info << "Thermal Diffusion Coefficient: D" << endl;
  Info << "Mass Heat Capacity: Cp" << endl;
  Info << "Chemical Enthalpy: dH" << endl;
  Info << "G-Yc,Yc: dCC" << endl;
  Info << "G-Yc,Z: dCZ" << endl;
  Info << "G-Z,Yc: dZC" << endl;
  Info << "G-Z,Z: dZZ" << endl;
  Info << "}\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

lookupFGM::~lookupFGM()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::lookupFGM::interpolateValue4D
(
    const List<scalar>& table,
    const scalar ZValue,
    const scalar CValue,
    const scalar ZvValue,
    const scalar CvValue,
    const List<scalar>& ZVec,
    const List<scalar>& CVec,
    const List<scalar>& ZvVec,
    const List<scalar>& CvVec,
    const int& nC,
    const int& nZv,
    const int& nCv
) const
{
    const scalar smallValue = 1e-12;
    const label nEntries = table.size();
    const label Nz = nEntries / (nC * nZv * nCv);
	
    if (nEntries != Nz * nC * nZv * nCv)
    {
        FatalErrorInFunction
            << "Table size mismatch: " << nEntries
            << " != " << Nz << "*" << nC << "*" << nZv << "*" << nCv
            << exit(FatalError);
    }
	
    // --- clamp queries ---
    scalar zq  = std::max(ZVec[0],  std::min(ZValue,  ZVec[Nz-1]));
    scalar cq  = std::max(CVec[0],  std::min(CValue,  CVec[nC-1]));
    scalar zvq = std::max(ZvVec[0], std::min(ZvValue, ZvVec[nZv-1]));
    scalar cvq = std::max(CvVec[0], std::min(CvValue, CvVec[nCv-1]));
	
    // --- find indices ---
    auto findIndex = [&](const List<scalar>& vec, scalar val, label n) -> label
    {
        for (label m = 1; m < n; ++m)
        {
            if (val <= vec[m]) return m - 1;
        }
        return n - 2;
    };
    label iz  = findIndex(ZVec,  zq,  Nz);
    label ic  = findIndex(CVec,  cq,  nC);
    label izv = findIndex(ZvVec, zvq, nZv);
    label icv = findIndex(CvVec, cvq, nCv);
	
    // --- interpolation weights ---
    scalar tz  = (zq  - ZVec[iz])      / std::max(ZVec[iz+1]   - ZVec[iz],   smallValue);
    scalar tc  = (cq  - CVec[ic])      / std::max(CVec[ic+1]   - CVec[ic],   smallValue);
    scalar tzv = (zvq - ZvVec[izv])    / std::max(ZvVec[izv+1] - ZvVec[izv], smallValue);
    scalar tcv = (cvq - CvVec[icv])    / std::max(CvVec[icv+1] - CvVec[icv], smallValue);

    // --- indexing lambda ---
    // Z (slowest) -> C -> Zv -> Cv (fastest)
    auto idx = [&](label iz_, label ic_, label izv_, label icv_) -> label
    {
        return iz_  * (nC * nZv * nCv)
             + ic_  * (nZv * nCv)
             + izv_ * (nCv)
             + icv_;
    };

    // --- fetch 16 corner values ---
    // f[dz][dc][dzv][dcv]
    scalar f[2][2][2][2];
    for (int dz  = 0; dz  < 2; ++dz)
    for (int dc  = 0; dc  < 2; ++dc)
    for (int dzv = 0; dzv < 2; ++dzv)
    for (int dcv = 0; dcv < 2; ++dcv)
    {
        f[dz][dc][dzv][dcv] =
            table[idx(iz+dz, ic+dc, izv+dzv, icv+dcv)];
    }

    // --- quadrilinear interpolation, step by step ---
    // Step 1: collapse Zv dimension
    scalar f1[2][2][2];
    for (int dz  = 0; dz  < 2; ++dz)
    for (int dc  = 0; dc  < 2; ++dc)
    for (int dcv = 0; dcv < 2; ++dcv)
    {
        f1[dz][dc][dcv] =
            f[dz][dc][0][dcv] +
            tzv * (f[dz][dc][1][dcv] - f[dz][dc][0][dcv]);
    }
	
    // Step 2: collapse Cv dimension
    scalar f2[2][2];
    for (int dz = 0; dz < 2; ++dz)
    for (int dc = 0; dc < 2; ++dc)
    {
        f2[dz][dc] =
            f1[dz][dc][0] +
            tcv * (f1[dz][dc][1] - f1[dz][dc][0]);
    }
	
    // Step 3: collapse C dimension
    scalar f3[2];
    for (int dz = 0; dz < 2; ++dz)
    {
        f3[dz] =
            f2[dz][0] +
            tc * (f2[dz][1] - f2[dz][0]);
    }
	
    // Step 4: collapse Z dimension
    scalar fFinal = f3[0] + tz * (f3[1] - f3[0]);
	
    return fFinal;
}

Foam::scalar Foam::lookupFGM::interpolateValue1D
(
    const List<scalar>& table,   
    const scalar zValue,         
    const List<scalar>& Zvec     
) const
{
	scalar smallValue = 1e-12;

    const label Nz = Zvec.size();

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

