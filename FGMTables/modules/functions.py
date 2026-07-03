# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:40:16 2026

@author: anton
"""


###############################################################
#
# FUNCTIONS 
#
###############################################################


#import :
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import matplotlib as mpl

mpl.rcParams['text.usetex'] = True   # enable LaTeX
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern']


def BilgerMixtureFraction(gas,f,mech):
    
    gas1 = ct.Solution(mech)
    air = "O2:0.21,N2:0.79"
    gas1.set_equivalence_ratio(phi=0, fuel="H2:1", oxidizer=air)
    gas1.TPX = 298.15, 101325, gas1.X
    
    
    W_H = gas.molecular_weights[gas.species_index("H")]
    W_O = gas.molecular_weights[gas.species_index("O")]
      
    Z_H = f.elemental_mass_fraction("H")
    Z_O = f.elemental_mass_fraction("O")

    Z_num = 1/(2*W_H)*(Z_H)-1/W_O*(Z_O-gas1.Y[gas.species_index('O2')])
    Z_den = 1/(2*W_H)-1/W_O*(-gas1.Y[gas.species_index('O2')])
    Z = Z_num/Z_den
    
    return Z
  


def premixedFlame(p,tin,PHI,mech,Soret):

    
    gas = ct.Solution(mech)

    air = "O2:0.21,N2:0.79"
    gas.set_equivalence_ratio(phi=PHI, fuel="H2:1", oxidizer=air)
    gas.TPX = tin, p, gas.X
    #Initial grid, chosen to be 0.02cm long :
    # - Refined grid at inlet and outlet, 6 points in x-direction :
    # initial_grid = 2*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/3
    initial_grid = np.linspace(0, 0.03, 100)
    #Create the free laminar premixed flame
    f = ct.FreeFlame(gas, initial_grid)
    max_grid_points = 10000
    f.max_time_step_count = 10000

    #First flame:
    #No energy for starters
    f.set_max_grid_points(f.flame, max_grid_points)
    #set inlet conditions  
    f.inlet.X = gas.X
    f.inlet.T = tin

    #First flame:
    #No energy for starters
    f.energy_enabled = False
    refine_grid = False
    #Tolerance properties
    tol_ss = [1.0e-5, 1.0e-9] # [rtol atol] for steady-state problem
    tol_ts = [1.0e-5, 1.0e-9] # [rtol atol] for time stepping
    # tol_ss = [1.0e-7, 1.0e-12] # [rtol atol] for steady-state problem
    # tol_ts = [1.0e-7, 1.0e-12] # [rtol atol] for time stepping
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    #Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(10, 10)
    #Set time steps whenever Newton convergence fails
    f.set_time_step(1.0e-5, [2, 5, 10, 20, 40, 50]) #s
    #Refinement criteria
    f.set_refine_criteria(ratio = 10.0, slope = 1, curve = 1)
    
    #Calculation
    loglevel = 0 # amount of diagnostic output (0 to 5)
    
    # gas.transport_model = 'unity-Lewis-number'
    # f.transport_model = 'unity-Lewis-number'
    gas.transport_model = 'mixture-averaged'
    f.transport_model = 'mixture-averaged'
    # gas.transport_model = 'multicomponent'
    # f.transport_model = 'multicomponent'
    
    f.solve(loglevel, refine_grid)
    refine_grid = True 
    f.energy_enabled = True
    f.soret_enabled = Soret
    #Refinement criteria when energy equation is enabled
    f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
    #Calculation and save of the results
    f.solve(loglevel, refine_grid)
    # f.save('ch4_adiabatic.yaml','energy','solution with the energy equation enabled')
    #See the sl to get an idea of whether or not you should continue
    points = f.flame.n_points
    #print('First Flame, mixture-averaged flamespeed = {0:7f} m/s'.format(f.velocity[0])) #m/s
    #print('First Flame, final T = {0:7f} K'.format(f.T[points-1])) #K
    #Second flame:
    #Energy equation enabled
    f.energy_enabled = True
    refine_grid = True
    f.set_refine_criteria(ratio = 2.0, slope = 0.01, curve = 0.001, prune = 1e-8)
    f.solve(loglevel, refine_grid)
    points = f.flame.n_points
    #print('Second Flame. laminar flamespeed = {0:7f} m/s'.format(f.velocity[0])) #m/s
    #print('Second Flame, final T = {0:7f} K'.format(f.T[points-1])) #K

    return f, gas


def phi2Z(phi,Y_fuel_u,Y_ox_u):
    
    s = 8
    Z = phi*Y_ox_u/(s+phi*Y_ox_u)
    
    return Z

def Z2phi(Z,Y_fuel_u,Y_ox_u):
    
    s = 8
    phi = s*Z/(Y_ox_u*(1-Z))
    
    return phi


def computeConsumptionSpeed(f, T_burnt, fuel, mech):
    """
    Compute the laminar consumption speed. Tho definitions are used, ie
    based on the integral of the heat release rate (Extended TFC model)
    and based on the 'fuel' species net consumption rate (suitable for Klarmann model)
    """
             
    rho_u = f.density[0]           
    
    #alternative based on fuel-species net production rate
    Y = 0
    Y2 = 0
    omega = 0
    for spec in fuel:
    	index = f.gas.species_index(fuel)
    	M = f.gas.molecular_weights[index]      # molecular weight [kg/kmol]
    	Y += f.Y[index][0]                       # mass fraction at fresh inlet
    	Y2 += f.Y[index][-1]                     # mass fraction at outlet   (added from cuenot tutorial)
    	omega += M*f.net_production_rates[index]   # net production rate [kg/m3 s] 
    
    Sc = (1/(rho_u*(Y-Y2)))*(-np.trapz(omega, f.grid))
    return Sc


def interpVar(var,phi,nPointsC,results):
    
    c_tab = results[phi]["c"]
    Var = results[phi][var]
    c_new = np.linspace(0, 1, nPointsC)
    var_interp = np.interp(c_new, c_tab, Var)
    
    return var_interp


def SavitzkyGolay(X,Y, direction, wind, order):
    
    window = wind        
    polyorder = order   
    
    dX = savgol_filter(X, window, polyorder, deriv=1, axis=direction)
    dY = savgol_filter(Y, window, polyorder, deriv=1, axis=direction)
    
    eps = 1e-12
    
    return dX / (dY + eps)


def interpExt(vecExt, vec, var):

    f = interp1d(vec, var, fill_value='extrapolate')
    return f(vecExt)


def interpExt2(vecExt, vec, var, varAir, varFuel):

    # Sort 
    idx = np.argsort(vec)
    Z_sorted = vec[idx]
    var_sorted = var[idx]

    # Remove duplicates by averaging values
    Z_unique, indices = np.unique(Z_sorted, return_inverse=True)
    var_unique = np.zeros_like(Z_unique, dtype=float)
    counts = np.zeros_like(Z_unique, dtype=int)

    for i, ind in enumerate(indices):
        var_unique[ind] += var_sorted[i]
        counts[ind] += 1

    var_unique /= counts 
    
    Z_unique = np.concatenate(([0], Z_unique, [1]))
    var_unique = np.concatenate(([varAir], var_unique, [varFuel]))
   
    # Interpolation
    f = interp1d(Z_unique, var_unique, fill_value='extrapolate')
    
    return f(vecExt)


def createTables(phis,var,results,nPointsC):
    
    Var_matrix = np.zeros((nPointsC,np.size(phis)))
    
    for i,phi in enumerate(phis):

        Var = results[phi][var]
        Var_matrix[:,i] = Var
    
    Var_matrix = np.array(Var_matrix)
    
    return Var_matrix  


def plot3DTable(C,Z,Table,SI,varTitle,xLab,yLab):
    
    fig = plt.figure(figsize=(6, 4), dpi=200, constrained_layout=True)
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(C, Z, Table, cmap='inferno', alpha=1)
    fig.colorbar(surf, location='top', shrink=0.5, aspect=20,
                  label=r'$'+varTitle+' '+SI+'$')

    # plt.title(r'$'+varTitle+' '+SI+'$')
    plt.xlabel(r'$'+xLab+'$')
    plt.ylabel(r'$'+yLab+'$')
    ax.set_box_aspect([2,2,2])  
    
def plotContourTable(C,Z,Table,SI,varTitle,xLab,yLab):
    
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize = (6,4), dpi = 400)
    contour = plt.contourf(C, Z, Table, 40, cmap='inferno')  
    plt.colorbar(contour, label=r'$'+varTitle+' '+SI+'$')
    # plt.title()
    plt.xlabel(r'$'+xLab+'$')
    plt.ylabel(r'$'+yLab+'$')
    plt.tight_layout()
    
    
def calculateLambdaYkC(Dk,Yk,W,C,species,speciesk,window,order):
    
    S = 0
    S1 = 0
    for sp in species:
        S += Dk[sp]*SavitzkyGolay(Yk[sp], C, 0, window , order)
        S1 += Dk[sp]*Yk[sp]
        
    dYkdC = SavitzkyGolay(Yk[speciesk], C, 0, window , order)
    dWdC = SavitzkyGolay(W, C, 0, window , order)
    LambdaYkC = Dk[speciesk]*dYkdC + (Dk[speciesk]*Yk[speciesk])/W*dWdC - Yk[speciesk]*S - Yk[speciesk]/W*dWdC*S1

    return LambdaYkC

def calculateLambdaYkZ(Dk,Yk,W,Z,species,speciesk,window,order):
    
    S = 0
    S1 = 0
    for sp in species:
        S += Dk[sp]*SavitzkyGolay(Yk[sp], Z, 1, window , order)
        S1 += Dk[sp]*Yk[sp]
        
    dYkdZ = SavitzkyGolay(Yk[speciesk], Z, 1, window , order)
    dWdZ = SavitzkyGolay(W, Z, 1, window , order)
    LambdaYkZ = Dk[speciesk]*dYkdZ + (Dk[speciesk]*Yk[speciesk]/W)*dWdZ - Yk[speciesk]*S - Yk[speciesk]/W*dWdZ*S1

    return LambdaYkZ

def calculateLambdaYkCSoret(Dk,DkSoret,T,rho,Yk,W,C,species,speciesk,window,order):
    
    S = 0
    S1 = 0
    S2 = 0
    for sp in species:
        S += Dk[sp]*SavitzkyGolay(Yk[sp], C, 0, window , order)
        S1 += Dk[sp]*Yk[sp]
        S2 += DkSoret[sp]/rho/T*SavitzkyGolay(T, C, 0, window , order)
        
    dYkdC = SavitzkyGolay(Yk[speciesk], C, 0, window , order)
    dWdC = SavitzkyGolay(W, C, 0, window , order)
    dTdC = SavitzkyGolay(T, C, 0, window , order)
    LambdaYkC = Dk[speciesk]*dYkdC + (Dk[speciesk]*Yk[speciesk])/W*dWdC - Yk[speciesk]*S - Yk[speciesk]/W*dWdC*S1 + DkSoret[speciesk]/rho/T*dTdC - Yk[speciesk]*S2

    return LambdaYkC

def calculateLambdaYkZSoret(Dk,DkSoret,T,rho,Yk,W,Z,species,speciesk,window,order):
    
    S = 0
    S1 = 0
    S2 = 0
    for sp in species:
        S += Dk[sp]*SavitzkyGolay(Yk[sp], Z, 1, window , order)
        S1 += Dk[sp]*Yk[sp]
        S2 += DkSoret[sp]/rho/T*SavitzkyGolay(T, Z, 1, window , order)
        
    dYkdZ = SavitzkyGolay(Yk[speciesk], Z, 1, window , order)
    dWdZ = SavitzkyGolay(W, Z, 1, window , order)
    dTdZ = SavitzkyGolay(T, Z, 1, window , order)
    LambdaYkZ = Dk[speciesk]*dYkdZ + (Dk[speciesk]*Yk[speciesk]/W)*dWdZ - Yk[speciesk]*S - Yk[speciesk]/W*dWdZ*S1 + DkSoret[speciesk]/rho/T*dTdZ - Yk[speciesk]*S

    return LambdaYkZ


def calculateLambdaYkh(Dk,Yk,W,h,species,speciesk,window,order):
    
    S = 0
    S1 = 0
    
    for sp in species:
        S += Dk[sp]*SavitzkyGolay(Yk[sp], h, 0, window , order)
        S1 += Dk[sp]*Yk[sp]
        
    dYkdh = SavitzkyGolay(Yk[speciesk], h, 0, window , order)
    dWdh = SavitzkyGolay(W, h, 0, window , order)
    LambdaYkh = Dk[speciesk]*dYkdh + (Dk[speciesk]*Yk[speciesk])/W*dWdh - Yk[speciesk]*S - Yk[speciesk]/W*dWdh*S1

    return LambdaYkh


def addBoundaries3D(Matrix, laminarTable, CTable, CTilde01, nZTilde, nCTilde, nCVar):
    
   Matrix01 = np.zeros((nZTilde, nCTilde, nCVar))
   Matrix01[:,1:-1,1:-1] = Matrix
   
   for i in range(nZTilde):
       
       Matrix01[i,0,:] = np.ones((nCVar))*laminarTable[0,i]
       Matrix01[i,-1,:] = np.ones((nCVar))*laminarTable[-1,i]
       Matrix01[i,:,0] = np.interp(CTilde01, CTable[:,0], laminarTable[:,i])
       Matrix01[i,1:-1,-1] = Matrix[i,:,-1]

   return Matrix01


def addBoundaries4D(Matrix, MatrixIntOverC01, laminarTable, CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar):
    
    Matrix01 = np.zeros((nZTilde, nZVar, nCTilde, nCVar))
    Matrix01[1:-1,1:-1,:,:] = Matrix
    Matrix01[:,0,:,:] = MatrixIntOverC01
    
    for i in range(nZVar):

        Matrix01[0,i,:,:] = MatrixIntOverC01[0,:,:]
        Matrix01[-1,i,:,:] = MatrixIntOverC01[-1,:,:]
     
    Matrix01[:,-1,:,:] = Matrix01[:,-2,:,:]

    return Matrix01

 
def writeFGMTable4D(
                        zLean, 
                        zRich, 
                        T0,
                        p,
                        rhoAir, 
                        kAir, 
                        CpAir, 
                        WAir, 
                        DAir, 
                        DNOAir,
                        rhoFuel, 
                        kFuel, 
                        CpFuel, 
                        WFuel, 
                        DFuel,
                        DNOFuel,
                        Z, 
                        C,
                        Zv,
                        Cv,
                        sL,
                        tau,
                        lF,
                        sourceH2,
                        sourceNO,
                        sourcePV,
                        sourcePVV,
                        T, 
                        rho,
                        D,
                        DNO,
                        k, 
                        Cp, 
                        W, 
                        HRR, 
                        dH,
                        YH, 
                        YH2, 
                        YH2O, 
                        YH2O2, 
                        YHO2, 
                        YO, 
                        YO2, 
                        YOH, 
                        YN2, 
                        YNO,
                        gYcYc, 
                        gYcZ, 
                        gZYc, 
                        gZZ,  
                        nZ,
                        nC,
                        nZVar,
                        nCVar,
                        mech,
                        filename="fgmProperties"
                        ):
    
    header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  9                                |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fgmProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

    ender = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    
    zLean = "zLean zLean [0 0 0 0 0 0 0] "+str(zLean)+"; // Lean mixture fraction limit"
    zRich = "zRich zRich [0 0 0 0 0 0 0] "+str(zRich)+"; // Rich mixture fraction limit"
    
    T0 = "T0 T0 [0 0 0 1 0 0 0] "+str(T0)+"; // Unburnt temperature"
    R0 = "R0 R0 [1 2 -2 -1 -1 0 0] 8.314462; // Universal gas constant"
    pRef = "pRef pRef [1 -1 -2 0 0 0 0] "+str(p)+"; // Thermodynamic reference pressure"
    
    sLAir = "sLAir sLAir [0 1 -1 0 0 0 0] 0; // Laminar flame speed value in pure air"
    tauAir = "tauAir tauAir [0 0 0 0 0 0 0] 0; // Heat release factor value in pure air"
    lFAir = "lFAir lFAir [0 1 0 0 0 0 0] 0; // Laminar flame thickness value in pure air"
    rhoAir = "rhoAir rhoAir [1 -3 0 0 0 0 0] "+str(rhoAir)+"; // Cp value in pure air"
    kAir = "kAir kAir [1 1 -3 -1 0 0 0] "+str(kAir)+"; // Thermal conductivity value in pure air"
    CpAir = "CpAir CpAir [0 2 -2 -1 0 0 0] "+str(CpAir)+"; // Cp value in pure air"
    dHAir = "dHAir dHAir [0 2 -2 0 0 0 0] 0; // Enthalpy of formation value in pure air"
    WAir = "WAir WAir [1 0 0 0 -1 0 0] "+str(WAir)+"; // Mixture molecular weight value in pure air"
    # HRRAir = "HRRAir HRRAir [1 -1 -3 0 0 0 0] 0; // Heat Release Rate value in pure air"
    DAir = "DAir DAir [0 2 -1 0 0 0 0] "+str(DAir)+"; // Mixture diffusion coefficient value in pure air"
    DNOAir = "DNOAir DNOAir [0 2 -1 0 0 0 0] "+str(DNOAir)+"; // Mixture diffusion coefficient value in pure air"
    YHAir = "YHAir YHAir [0 0 0 0 0 0 0] 0; // H mass fraction value in pure air"
    YH2Air = "YH2Air YH2Air [0 0 0 0 0 0 0] 0; // H2 mass fraction value in pure air"
    YH2OAir = "YH2OAir YH2OAir [0 0 0 0 0 0 0] 0; // H2O mass fraction value in pure air"
    YH2O2Air = "YH2O2Air YH2O2Air [0 0 0 0 0 0 0] 0; // H2O2 mass fraction value in pure air"
    YHO2Air = "YHO2Air YHO2Air [0 0 0 0 0 0 0] 0; // HO2 mass fraction value in pure air"
    YOAir = "YOAir YOAir [0 0 0 0 0 0 0] 0; // O mass fraction value in pure air"
    YO2Air = "YO2Air YO2Air [0 0 0 0 0 0 0] 0.23290921795842306; // O2 mass fraction value in pure air"
    YOHAir = "YOHAir YOHAir [0 0 0 0 0 0 0] 0; // OH mass fraction value in pure air"
    YN2Air = "YN2Air YN2Air [0 0 0 0 0 0 0] 0.767090782041577; // N2 mass fraction value in pure air"
    # YOHEXAir = "YOHEXAir YOHEXAir [0 0 0 0 0 0 0] 0; // OH* mass fraction value in pure air"
    dCCAir = "dCCAir dCCAir [1 -1 -1 0 0 0 0] 0; // Gamma_{Yc,Yc}*rho value in pure Air"
    dCZAir = "dCZAir dCZAir [1 -1 -1 0 0 0 0] 0; // Gamma_{Yc,Z}*rho value in pure Air"
    dZCAir = "dZCAir dZCAir [1 -1 -1 0 0 0 0] 0; // Gamma_{Z,Yc}*rho value in pure Air"
    dZZAir =  "dZZAir dZZAir [1 -1 -1 0 0 0 0] 0; // Gamma_{Z,Z}*rho value in pure Air"
    
    sLFuel = "sLFuel sLFuel [0 1 -1 0 0 0 0] 0; // Laminar flame speed value in pure Fuel"
    tauFuel = "tauFuel tauFuel [0 0 0 0 0 0 0] 0; // Heat release factor value in pure Fuel"
    lFFuel = "lFFuel lFFuel[0 1 0 0 0 0 0] 0; // Laminar flame thickness value in pure Fuel"
    rhoFuel = "rhoFuel rhoFuel [1 -3 0 0 0 0 0] "+str(rhoFuel)+"; // Cp value in pure Fuel"
    kFuel = "kFuel kFuel [1 1 -3 -1 0 0 0] "+str(kFuel)+"; // Thermal conductivity value in pure Fuel"
    CpFuel = "CpFuel CpFuel [0 2 -2 -1 0 0 0] "+str(CpFuel)+"; // Cp value in pure Fuel"
    dHFuel = "dHFuel dHFuel [0 2 -2 0 0 0 0] 0; // Enthalpy of formation value in pure Fuel"
    WFuel = "WFuel WFuel [1 0 0 0 -1 0 0] "+str(WFuel)+"; // Mixture molecular weight value in pure Fuel"
    # HRRFuel = "HRRFuel HRRFuel [1 -1 -3 0 0 0 0] 0; // Heat Release Rate value in pure Fuel"
    DFuel = "DFuel DFuel [0 2 -1 0 0 0 0] "+str(DFuel)+"; // Mixture diffusion coefficient value in pure Fuel"
    DNOFuel = "DNOFuel DNOFuel [0 2 -1 0 0 0 0] "+str(DNOFuel)+"; // Mixture diffusion coefficient value in pure Fuel"

    YHFuel = "YHFuel YHFuel [0 0 0 0 0 0 0] 0; // H mass fraction value in pure Fuel"
    YH2Fuel = "YH2Fuel YH2Fuel [0 0 0 0 0 0 0] 1; // H2 mass fraction value in pure Fuel"
    YH2OFuel = "YH2OFuel YH2OFuel  [0 0 0 0 0 0 0] 0; // H2O mass fraction value in pure Fuel"
    YH2O2Fuel = "YH2O2Fuel YH2O2Fuel [0 0 0 0 0 0 0] 0; // H2O2 mass fraction value in pure Fuel"
    YHO2Fuel = "YHO2Fuel YHO2Fuel [0 0 0 0 0 0 0] 0; // HO2 mass fraction value in pure Fuel"
    YOFuel = "YOFuel YOFuel [0 0 0 0 0 0 0] 0; // O mass fraction value in pure Fuel"
    YO2Fuel = "YO2Fuel YO2Fuel [0 0 0 0 0 0 0] 0; // O2 mass fraction value in pure Fuel"
    YOHFuel = "YOHFuel YOHFuel [0 0 0 0 0 0 0] 0; // OH mass fraction value in pure Fuel"
    YN2Fuel = "YN2Fuel YN2Fuel [0 0 0 0 0 0 0] 0; // N2 mass fraction value in pure Fuel"
    # YOHEXFuel = "YOHEXFuel YOHEXFuel [0 0 0 0 0 0 0] 0; // OH* mass fraction value in pure Fuel"
    dCCFuel = "dCCFuel dCCFuel [1 -1 -1 0 0 0 0] 0; // Gamma_{Yc,Yc}*rho value in pure Fuel"
    dCZFuel = "dCZFuel dCZFuel [1 -1 -1 0 0 0 0] 0; // Gamma_{Yc,Z}*rho value in pure Fuel"
    dZCFuel = "dZCFuel dZCFuel [1 -1 -1 0 0 0 0] 0; // Gamma_{Z,Yc}*rho value in pure Fuel"
    dZZFuel =  "dZZFuel dZZFuel [1 -1 -1 0 0 0 0] 0; // Gamma_{Z,Z}*rho value in pure Fuel"
    
    nPointsZ = "nZ "+str(nZ)+";"
    nPointsC = "nC "+str(nC)+";"
    nPointsZVar = "nZv "+str(nZVar)+";"
    nPointsCVar = "nCv "+str(nCVar)+";"
    
    with open(filename, "w") as f:
        f.write(
            header + "\n\n" + 
            "//Reaction Mechanism: " + mech + "\n\n" +
            "//Table Discretization" + "\n\n" +
            nPointsZ + "\n" + 
            nPointsC + "\n" + 
            nPointsZVar + "\n" + 
            nPointsCVar + "\n\n" + 
            "//Flammability Limits" + "\n\n" +
            zLean + "\n" +
            zRich + "\n\n" +
            "//Thermodynamic Reference Values" + "\n\n" +
            T0 + "\n" +
            R0 + "\n" +
            pRef + "\n\n" +
            "//Values in Pure Air" + "\n\n" +
            sLAir + "\n" +
            tauAir + "\n" +
            lFAir + "\n" +
            rhoAir + "\n" +
            kAir + "\n" +
            CpAir + "\n" +
            dHAir + "\n" +
            WAir + "\n" +
            DAir + "\n" +
            DNOAir + "\n" +
            YHAir + "\n" +
            YH2Air + "\n" +
            YH2OAir + "\n" +
            YH2O2Air + "\n" +
            YHO2Air + "\n" +
            YOAir + "\n" +
            YO2Air + "\n" +
            YOHAir + "\n" +
            YN2Air + "\n" +
            dCCAir + "\n" +
            dCZAir + "\n" +
            dZCAir + "\n" +
            dZZAir + "\n\n" +
            "//Values in Pure Fuel" + "\n\n" +
            sLFuel + "\n" +
            tauFuel + "\n" +
            lFFuel + "\n" +
            rhoFuel + "\n" +
            kFuel + "\n" +
            CpFuel + "\n" +
            dHFuel + "\n" +
            WFuel + "\n" +
            DFuel + "\n" +
            DNOFuel + "\n" +
            YHFuel + "\n" +
            YH2Fuel + "\n" +
            YH2OFuel + "\n" +
            YH2O2Fuel + "\n" +
            YHO2Fuel + "\n" +
            YOFuel + "\n" +
            YO2Fuel + "\n" +
            YOHFuel + "\n" +
            YN2Fuel + "\n" +
            dCCFuel + "\n" +
            dCZFuel + "\n" +
            dZCFuel + "\n" +
            dZZFuel + "\n\n"
            )
        
        f.write("//Tabulated Values Z - C - Cv - Zv\n\n")

        write_vector(f, 'Z', Z)
        write_vector(f, 'C', C)
        write_vector(f, 'Zv', Zv)
        write_vector(f, 'Cv', Cv)
        write_vector(f, 'sL', sL)
        write_vector(f, 'tau', tau)
        write_vector(f, 'lF', lF)
        write_vector(f, 'omegaYc', sourcePV)
        write_vector(f, 'omegaYcV', sourcePVV)
        write_vector(f, 'omegaH2', sourceH2)
        write_vector(f, 'omegaNO', sourceNO)
        write_vector(f, 'T', T)
        write_vector(f, 'Cp', Cp)
        write_vector(f, 'dH', dH)
        write_vector(f, 'D', D)
        write_vector(f, 'DNO', DNO)
        write_vector(f, 'rho', rho)
        write_vector(f, 'k', k)
        write_vector(f, 'W', W)
        write_vector(f, 'HRR', HRR)
        write_vector(f, 'Yk_H', YH)
        write_vector(f, 'Yk_H2', YH2)
        write_vector(f, 'Yk_H2O', YH2O)
        write_vector(f, 'Yk_H2O2', YH2O2)
        write_vector(f, 'Yk_HO2', YHO2)
        write_vector(f, 'Yk_O', YO)
        write_vector(f, 'Yk_O2', YO2)
        write_vector(f, 'Yk_OH', YOH)
        write_vector(f, 'Yk_N2', YN2)
        write_vector(f, 'Yk_NO', YNO)
        write_vector(f, 'dYcYc', gYcYc)
        write_vector(f, 'dYcZ', gYcZ)
        write_vector(f, 'dZYc', gZYc)
        write_vector(f, 'dZZ', gZZ)
        
        f.write(ender)

    print(f"File '{filename}' successfully created.")
    
def write_vector(f, name, values):
    f.write(f"{name}\n(\n")
    for v in values:
        f.write(f"{v}\n")
    f.write(");\n\n")
    
    
    
def scaleYc(Z,YcMin,YcMax, filename="scalingYcTable"):
    header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  9                                |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      scalingYcTable;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

    ender = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    

    def format_vector(name, values):
        formatted_values = "\n".join(f"{v}" for v in values)
        return f"{name}\n(\n{formatted_values}\n);"

    body = (
            f"{format_vector('ZScal', Z)}\n\n"
            f"{format_vector('YcMin', YcMin)}\n\n"
            f"{format_vector('YcMax', YcMax)}"
           )
    
    with open(filename, "w") as f:
        f.write(
                header + "\n" + 
                body + "\n\n" + 
                ender
            )

    print(f"File '{filename}' successfully created.")