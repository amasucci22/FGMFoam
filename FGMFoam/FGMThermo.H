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

Class
    Foam::heThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    heThermo.C

\*---------------------------------------------------------------------------*/

#ifndef FGMThermo_H
#define FGMThermo_H

#include "psiThermo.H"
#include "basicMixture.H"
#include "volFields.H"
#include "uniformGeometricFields.H"

#include "lookupFGM.H"
#include "scaleYc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class heThermo Declaration
\*---------------------------------------------------------------------------*/

class FGMThermo
:
    public psiThermo::composite
{
protected:

    // Protected data
	
		//- Look up table and routines
        const lookupFGM& FGMTable_;
	   
	    //- Look up table for Yc scaling
        const scaleYc& scaleYcTable_;
		
		//- Transport equations
		volScalarField& Yc_;
        volScalarField& Z_;
		volScalarField& YcV_;
        volScalarField& ZV_;

        //- Energy field
        volScalarField he_;

        //- Heat capacity at constant pressure field [J/kg/K]
        volScalarField Cp_;

        //- Heat capacity at constant volume field [J/kg/K]
        volScalarField Cv_;
		
		//- Enthalpy of formation [m^2/s^2]
		volScalarField dH_;
		
		//- Temperature [K]
		volScalarField Temp_;
	   
		//- YcMin [-]
		volScalarField YcMin_;
	   
		//- YcMax [-]
		volScalarField YcMax_;
	   
		//- YcMin [-]
		volScalarField ZScal_;
	   
		//- Scaled Progress variable [-]
		volScalarField C_;
 
		//- Rate reaction of progress variable [kg/(m^3/s)]
		volScalarField sourceYc_;
		
		//- Rate reaction of progress variable [kg/(m^3/s)]
		volScalarField sourceYcV_;
	   
		//- Rate reaction of H2 [kg/(m^3/s)]
		volScalarField sourceH2_;
		
		//- Rate reaction of NO [kg/(m^3/s)]
		volScalarField sourceNO_;					   
	   
		//- Gamma Yc-Yc [kg/(m s)]
		volScalarField dCC_;
	   
		//- Gamma Yc-Z [kg/(m s)]
		volScalarField dCZ_;
	   
		//- Gamma Z-Yc [kg/(m s)]
		volScalarField dZC_;
	   
		//- Gamma Z-Z [kg/(m s)]
		volScalarField dZZ_;
	   
		//- Thermal Conductivity [W/(mK)]
		volScalarField k_;
	   
		//- Mixture Diffusion Coefficient [m^2/s]
		volScalarField D_;
		
		//- NO Diffusion Coefficient [m^2/s]
		volScalarField DNO_;						   
	   
		//- Molecular Weight [kg/mol]
		volScalarField W_;
	   
		//- H Mass Fraction [-]
		volScalarField H_;
	   
		//- H2 Mass Fraction [-]
		volScalarField H2_;  
	   
		//- H2O Mass Fraction [-]
		volScalarField H2O_;
	   
		//- HO2 Mass Fraction [-]
		volScalarField HO2_;
	   
		//- H2O2 Mass Fraction [-]
		volScalarField H2O2_;
	   
		//- O Mass Fraction [-]
		volScalarField O_;
	   
		//- O2 Mass Fraction [-]
		volScalarField O2_;  
	
		//- OH Mass Fraction [-]
		volScalarField OH_;
	   
		//- N2 Mass Fraction [-]
		volScalarField N2_;
		
		//- NO Mass Fraction [-]
		volScalarField NO_;				  
	   
		//- Heat Release Rate [W/m^3]
		volScalarField HRR_;
		
		//- Laminar Flame Speed [m/s]
		volScalarField sL_;
		
		//- Laminar Flame Thickness [m]
		volScalarField lF_;
		
		//- Laminar Heat Release Factor [-]
		volScalarField tau_;
		
		//- Variance of C [-]
		volScalarField CV_;
		

    // Protected Member Functions

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);


public:

    // Constructors

        //- Construct from mesh
        FGMThermo
        (
		    const word&,
            const fvMesh&,
			const lookupFGM&,
			const scaleYc&,
			volScalarField&,
			volScalarField&,
			volScalarField&,
			volScalarField&,
			const word& phaseName=word::null
        );



    //- Destructor
    virtual ~FGMThermo();


    // Member Functions


        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return "FGMThermo";
        }
		
		//- Calculate the thermo variables
        void correct();

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return false;
        }


        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                return he_;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                return he_;
            }

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const
            {
                return Cp_;
            }
			

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const
            {
                return Cv_;
            }
			
			//- SGS Karlovitz number
			tmp<volScalarField> Ka
			(
				const volScalarField& k,
				const volScalarField& delta
			) const;

			//- SGS Damkholer number
			tmp<volScalarField> Da
			(
				const volScalarField& k,
				const volScalarField& delta
			) const;
			
			tmp<volScalarField> C2
			(
				const volScalarField& k,
				const volScalarField& T,
				const volScalarField& delta
			);
			
			//- central term of the SGS dissipation for the mix. frac. variance
            tmp<volScalarField> zvSGSdiss
			(
				const volScalarField& muT,
				const volScalarField& delta
			);
			
			tmp<volScalarField> filter(const tmp<volScalarField>& phi);
			


        // Fields derived from thermodynamic state variables

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy [J/kg/K]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Absolute enthalpy [J/kg/K]
            virtual tmp<volScalarField> ha() const;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Absolute enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const;
			
			//- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> Temp() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> Temp(const label patchi) const;
			


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappa(const label patchi) const;

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;
			
			tmp<volScalarField> delta() const;
			
			tmp<volScalarField> rho() const
            {
			if (incompressible())
			{
				  const dimensionedScalar pRef(FGMTable_.lookup("pRef"));
				  return pRef*psi();
			}
			else
			{
		      return p_*psi();
	        }
            }

            //Additional functions from version 2.2 on

            tmp<scalarField> rho(const label patchi) const
            {

                if (incompressible())
			{
				const dimensionedScalar pRef(FGMTable_.lookup("pRef"));
				return pRef.value()*psi_.boundaryField()[patchi];
			}
			else
			{
				 return p_.boundaryField()[patchi]*psi_.boundaryField()[patchi];
			}
				}
			
	
		
			//- Minimum Yc [-]
			virtual tmp<volScalarField> YcMin() const
			{
			  return YcMin_;
			}
		
			// Return non-const access to the local Minimum Yc [-]
			virtual volScalarField& YcMin()
			{
			  return YcMin_;
			}
			
			//- Maximum Yc [-]
			virtual tmp<volScalarField> YcMax() const
			{
			  return YcMax_;
			}
			
			// Return non-const access to the local Maximum Yc [-]
			virtual volScalarField& YcMax()
			{
			  return YcMax_;
			}
			
			//- Z scaling [-]
			virtual tmp<volScalarField> ZScal() const
			{
			  return ZScal_;
			}
			
			// Return non-const access to the local Z scaling [-]]
			virtual volScalarField& ZScal()
			{
			  return ZScal_;
			}
			
			//- Scaled Progress Variable
			virtual tmp<volScalarField> C() const
			{
			  return C_;
			}
			
			// Return non-const access to the local Scaled Progress Variable
			virtual volScalarField& C()
			{
			  return C_;
			}

			
			//- Dynamic viscosity [kg/m/s]
			tmp<volScalarField> mu() const;
			
		  
			//- Rate reaction of progress variable [kg/(m^3*s)]
			virtual tmp<volScalarField> sourceYc() const
			{
			  return sourceYc_;
			}

			//- Return non-const access to the local rate reaction of
			// progress variable field [kg/(m^3*s)]
			virtual volScalarField& sourceYc()
			{
			  return sourceYc_;
			}
			
			//- Rate reaction of progress variable variance  [kg/(m^3*s)]
			virtual tmp<volScalarField> sourceYcV() const
			{
			  return sourceYcV_;
			}

			//- Return non-const access to the local rate reaction of
			// progress variable field variance [kg/(m^3*s)]
			virtual volScalarField& sourceYcV()
			{
			  return sourceYcV_;
			}
			
			//- H2 source Term [kg/(m^3*s)]
			virtual tmp<volScalarField> sourceH2() const
			{
			  return sourceH2_;
			}
			
			// Return non-const access to the local H2 source Term [kg/(m^3*s)]
			virtual volScalarField& sourceH2()
			{
			  return sourceH2_;
			}
			
			//- NO source Term [kg/(m^3*s)]
			virtual tmp<volScalarField> sourceNO() const
			{
			  return sourceNO_;
			}
			
			// Return non-const access to the local NO source Term [kg/(m^3*s)]
			virtual volScalarField& sourceNO()
			{
			  return sourceNO_;
			}			   
			//- Gamma Yc-Yc
			virtual tmp<volScalarField> dCC() const
			{
			  return dCC_;
			}
			
			//- Return non-const access to Gamma Yc-Yc
			virtual volScalarField& dCC()
			{
			  return dCC_;
			}

			//- Gamma Yc-Z
			virtual tmp<volScalarField> dCZ() const
			{
			  return dCZ_;
			}
			
			//- Return non-const access to Gamma Yc-Z
			virtual volScalarField& dCZ()
			{
			  return dCZ_;
			}

			//- Gamma Z-Yc
			virtual tmp<volScalarField> dZC() const
			{
			  return dZC_;
			}
			
			//- Return non-const access to Gamma Z-Yc
			virtual volScalarField& dZC()
			{
			  return dZC_;
			}
			
			//- Gamma Z-Z
			virtual tmp<volScalarField> dZZ() const
			{
			  return dZZ_;
			}
			
			//- Return non-const access to Gamma Z-Z
			virtual volScalarField& dZZ()
			{
			  return dZZ_;
			}
			
			//- Thermal Conductivity [W/(mK)]
			virtual tmp<volScalarField> k() const
			{
			  return k_;
			}
			
			//- Return non-const access to Thermal Conductivity [W/(mK)]
			virtual volScalarField& k()
			{
			  return k_;
			}

			//- Mass diffusivity [m^2/s]
			virtual tmp<volScalarField> D() const
			{
			  return D_;
			}

			//- Return non-const access to the local mass diffusivity field [m^2/s]
			virtual volScalarField& D()
			{
			  return D_;
			}
			
			//- NO Mass diffusivity [m^2/s]
			virtual tmp<volScalarField> DNO() const
			{
			  return DNO_;
			}

			//- Return non-const access to the local NO mass diffusivity field [m^2/s]
			virtual volScalarField& DNO()
			{
			  return DNO_;
			}				   
			
			//- H Mass Fraction [-]
			virtual tmp<volScalarField> H() const
			{
			  return H_;
			}
			
			// Return non-const access to the local H Mass Fraction [-]
			virtual volScalarField& H()
			{
			  return H_;
			}		

			//- H2 Mass Fraction [-]
			virtual tmp<volScalarField> H2() const
			{
			  return H2_;
			}
			
			// Return non-const access to the local H2 Mass Fraction [-]
			virtual volScalarField& H2()
			{
			  return H2_;
			}
			
			//- H2O Mass Fraction [-]
			virtual tmp<volScalarField> H2O() const
			{
			  return H2O_;
			}
			
			// Return non-const access to the local H2O Mass Fraction [-]
			virtual volScalarField& H2O()
			{
			  return H2O_;
			}
			
			//- HO2 Mass Fraction [-]
			virtual tmp<volScalarField> HO2() const
			{
			  return HO2_;
			}
			
			// Return non-const access to the local HO2 Mass Fraction [-]
			virtual volScalarField& HO2()
			{
			  return HO2_;
			}
			
			//- H2O2 Mass Fraction [-]
			virtual tmp<volScalarField> H2O2() const
			{
			  return H2O2_;
			}
			
			// Return non-const access to the local H2O2 Mass Fraction [-]
			virtual volScalarField& H2O2()
			{
			  return H2O2_;
			}
			
			//- O Mass Fraction [-]
			virtual tmp<volScalarField> O() const
			{
			  return O_;
			}
			
			// Return non-const access to the local O Mass Fraction [-]
			virtual volScalarField& O()
			{
			  return O_;
			}	
			
			//- O2 Mass Fraction [-]
			virtual tmp<volScalarField> O2() const
			{
			  return O2_;
			}
			
			// Return non-const access to the local O2 Mass Fraction [-]
			virtual volScalarField& O2()
			{
			  return O2_;
			}	
			
			//- OH Mass Fraction [-]
			virtual tmp<volScalarField> OH() const
			{
			  return OH_;
			}
			
			// Return non-const access to the local OH Mass Fraction [-]
			virtual volScalarField& OH()
			{
			  return OH_;
			}
			
			//- N2 Mass Fraction [-]
			virtual tmp<volScalarField> N2() const
			{
			  return N2_;
			}
			
			// Return non-const access to the local B2 Mass Fraction [-]
			virtual volScalarField& N2()
			{
			  return N2_;
			}
			//- NO Mass Fraction [-]
			virtual tmp<volScalarField> NO() const
			{
			  return NO_;
			}
			
			// Return non-const access to the local B2 Mass Fraction [-]
			virtual volScalarField& NO()
			{
			  return NO_;
			}			   

			//- Heat Release Rate [W/m^3]
			virtual tmp<volScalarField> HRR() const
			{
			  return HRR_;
			}
			
			// Return non-const access to the local Heat Release Rate [W/m^3]
			virtual volScalarField& HRR()
			{
			  return HRR_;
			}
			
			//- Heat Release Rate [W/m^3]
			virtual tmp<volScalarField> dH() const
			{
			  return dH_;
			}
			
			// Return non-const access to the local Heat Release Rate [W/m^3]
			virtual volScalarField& dH()
			{
			  return dH_;
			}
			
			//- Scaled Progress Variable
			virtual tmp<volScalarField> CV() const
			{
			  return CV_;
			}
			
			// Return non-const access to the local Scaled Progress Variable
			virtual volScalarField& CV()
			{
			  return CV_;
			}


			//- Read thermophysical properties dictionary
			virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
