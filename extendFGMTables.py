# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:40:16 2026

@author: anton
"""


##########################################################################
#
# EXTEND FGM TABLES OUT OF THE FLAMABILITY LIMITS Z \in [0,1]
#
##########################################################################


#import :
import cantera as ct
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mpl
from scipy.interpolate import interp1d
import modules.functions as m

mpl.rcParams['text.usetex'] = True   # enable LaTeX
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern']

# Path to reaction mechanism
ct.add_directory("C:/Program Files/Cantera/data")

# Do you want to save the Tables?
save        = True
plotTables  = False
Soret       = False

# Import reaction mechanism
nameMech = 'SanDiegoNO'
mech = nameMech+'.yaml'

if Soret:
    # Path where tables are stored
    path = "./" + nameMech + "Soret/"
    pathOut = "./" + nameMech + "SoretExtended/"
    print("Soret activated!")
else:
    # Path where tables are stored
    path = "./" + nameMech + "/"
    pathOut = "./" + nameMech + "Extended/"
    print("Soret disactivated!")
        

#Import the scaling function for Yc,b(Z) and Yc,u(Z)
YScal = pd.read_csv(path + "scaleYc" + nameMech + ".csv", header=None, skiprows=1).to_numpy()

#
gas = ct.Solution(mech)
species = gas.species_names


# Load Tables
DkTable = {}
DkSoretTable = {}
YkTable = {}
omegakTable = {}
XkTable = {}
hkTable = {}

CTable = pd.read_csv(path + "CTable.csv",header=None, skiprows=1).to_numpy()
ZTable = pd.read_csv(path + "ZTable.csv",header=None, skiprows=1).to_numpy()
TTable = pd.read_csv(path + "TTable.csv",header=None, skiprows=1).to_numpy()
sLVec = pd.read_csv(path + "sLVec.csv",header=None, skiprows=1).to_numpy()
tauVec = pd.read_csv(path + "tauVec.csv",header=None, skiprows=1).to_numpy()
lFVec = pd.read_csv(path + "lFVec.csv",header=None, skiprows=1).to_numpy()
dHTable = pd.read_csv(path + "dHTable.csv",header=None, skiprows=1).to_numpy()
DTable = pd.read_csv(path + "DTable.csv",header=None, skiprows=1).to_numpy()
rhoTable = pd.read_csv(path + "rhoTable.csv",header=None, skiprows=1).to_numpy()
kTable = pd.read_csv(path + "kTable.csv",header=None, skiprows=1).to_numpy()
CpTable = pd.read_csv(path + "CpTable.csv",header=None, skiprows=1).to_numpy()
MMWTable = pd.read_csv(path + "MMWTable.csv",header=None, skiprows=1).to_numpy()
HRRTable = pd.read_csv(path + "HRRTable.csv",header=None, skiprows=1).to_numpy()
hTable = pd.read_csv(path + "hTable.csv",header=None, skiprows=1).to_numpy()
sourceYcVTable = pd.read_csv(path + "sourceYcVTable.csv",header=None, skiprows=1).to_numpy()

for sp in species:
    omegakTable[sp] = pd.read_csv(path + f"omega_{sp}Table.csv",header=None, skiprows=1).to_numpy()
    DkTable[sp] = pd.read_csv(path + f"D_{sp}Table.csv",header=None, skiprows=1).to_numpy()
    DkSoretTable[sp] = pd.read_csv(path + f"DSoret_{sp}Table.csv",header=None, skiprows=1).to_numpy()
    YkTable[sp] = pd.read_csv(path + f"Y_{sp}Table.csv",header=None, skiprows=1).to_numpy()
    XkTable[sp] = pd.read_csv(path + f"X_{sp}Table.csv",header=None, skiprows=1).to_numpy()
    hkTable[sp] = pd.read_csv(path + f"h_{sp}Table.csv",header=None, skiprows=1).to_numpy()


ZVec = ZTable[0,:]
CVec = CTable[:,0]

zLean = ZVec[0]
zRich = ZVec[-1]


nPointsC = len(CVec)
nPointsZ = len(ZVec)

YcMax = np.zeros(len(ZVec))
YcMin = np.zeros(len(ZVec))


# Set Thermodynamic parameters for the flamelets
p = ct.one_atm  #pressure
tin = 298.15 # unburt gas temperature

# Create Air and Fuel properties 
# Air
gasAir = ct.Solution(mech)
air = "O2:0.21,N2:0.79"
gasAir.set_equivalence_ratio(phi=0, fuel="H2:1", oxidizer=air)
gasAir.TPX = tin, p, gasAir.X
YO2Air = gasAir.Y[gasAir.species_index('O2')]

Air = {
       "CpAir"      : gasAir.cp_mass,
       "WAir"       : gasAir.mean_molecular_weight*1e-3,
       "DAir"       : gasAir.thermal_conductivity/(gasAir.density*gasAir.cp_mass),
       "rhoAir"     : gasAir.density,
       "kAir"       : gasAir.thermal_conductivity,
       "HRRAir"     : 0
       }

for i, sp in enumerate(species):
    Air[f"D_{sp}Air"] = gasAir.mix_diff_coeffs[gasAir.species_index(sp)]
    Air[f"Y_{sp}Air"] = gasAir.Y[gasAir.species_index(sp)]
    Air[f"X_{sp}Air"] = gasAir.X[gasAir.species_index(sp)]

# Fuel
gasFuel = ct.Solution(mech)
gasFuel.TPX = tin, p, "H2:1"

Fuel = {
       "CpFuel"      : gasFuel.cp_mass,
       "WFuel"       : gasFuel.mean_molecular_weight*1e-3,
       "DFuel"       : gasFuel.thermal_conductivity/(gasFuel.density*gasFuel.cp_mass),
       "rhoFuel"     : gasFuel.density,
       "kFuel"       : gasFuel.thermal_conductivity,
       "HRRFuel"     : 0
       }

for i, sp in enumerate(species):
    Fuel[f"D_{sp}Fuel"] = gasFuel.mix_diff_coeffs[gasFuel.species_index(sp)]
    Fuel[f"Y_{sp}Fuel"] = gasFuel.Y[gasFuel.species_index(sp)]
    Fuel[f"X_{sp}Fuel"] = gasFuel.X[gasFuel.species_index(sp)]
   



# Create Z matrix where Z is constant along the rows and varies along the columns
ZTableExt1 = np.zeros((nPointsC,nPointsZ))
CTableExt = CTable

for i in range(nPointsZ):
    ZTableExt1[:,i] = np.ones(nPointsC)*ZTable[0,i]
    

# Interoplate over Z
for i in range(nPointsC):
    TTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], TTable[i,:])
    CpTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], CpTable[i,:])
    DTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], DTable[i,:])
    kTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], kTable[i,:])
    rhoTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], rhoTable[i,:])
    dHTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], dHTable[i,:])
    hTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], hTable[i,:])
    MMWTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], MMWTable[i,:])
    HRRTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], HRRTable[i,:])
    sourceYcVTable[i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], sourceYcVTable[i,:])
    
    for sp in species:
        YkTable[sp][i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], YkTable[sp][i,:])
        hkTable[sp][i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], hkTable[sp][i,:])
        XkTable[sp][i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], XkTable[sp][i,:])
        DkTable[sp][i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], DkTable[sp][i,:])
        DkSoretTable[sp][i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], DkSoretTable[sp][i,:])
        omegakTable[sp][i,:] = m.interpExt(ZTableExt1[i,:], ZTable[i,:], omegakTable[sp][i,:])
        
        
windowC = 40
windowZ = 40
order = 1



# Gamma_YcYc
gammaYcYc1 = np.zeros((nPointsC, nPointsZ))

if Soret:
    LYcC = m.calculateLambdaYkC(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, CTableExt, species, 'H2O', windowC, order)
else:
    LYcC = m.calculateLambdaYkC(DkTable, YkTable, MMWTable, CTableExt, species, 'H2O', windowC, order)

for i in range(nPointsZ):
    gammaYcYc1[:, i] = LYcC[:, i] / YkTable['H2O'][:, i][-1] - DTable[:, i]
gammaYcYc = (gammaYcYc1 + DTable) * rhoTable

# Gamma_YcZ
gammaYcZ  = np.zeros((nPointsC, nPointsZ))
dH2OdZ    = m.SavitzkyGolay(YkTable['H2O'], ZTableExt1, 1, windowC, order)

if Soret:
    LYcZ      = m.calculateLambdaYkZ(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, ZTableExt1, species, 'H2O', windowC, order)
else:
    LYcZ = m.calculateLambdaYkZ(DkTable, YkTable, MMWTable, ZTableExt1, species, 'H2O', windowC, order)

for i in range(nPointsC):
    dcdz_YcConst  = -1 / YkTable['H2O'][-1, :] * (CTableExt[i, :] * dH2OdZ[-1, :])
    gammaYcZ[i,:] = LYcC[i, :] * dcdz_YcConst + LYcZ[i, :]
gammaYcZ = gammaYcZ * rhoTable

# Gamma_ZYc
gammaZYc     = np.zeros((nPointsC, nPointsZ))

if Soret:
    LambdakTable = {sp: m.calculateLambdaYkC(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, CTableExt, species, sp, windowZ, order) for sp in species}
else:
    LambdakTable = {sp: m.calculateLambdaYkC(DkTable, YkTable, MMWTable, CTableExt, species, sp, windowZ, order) for sp in species}

gas = ct.Solution(mech)

listSpecies = ['H', 'O']
LZpCTable = {}
for k in listSpecies:
    isp = gas.element_index(k)
    S = sum(
        gas.n_atoms(gas.species_index(sp), isp) / gas.molecular_weights[gas.species_index(sp)] * LambdakTable[sp]
        for sp in species
    )
    LZpCTable[k] = gas.molecular_weights[gas.species_index(k)] * S

KZ  = 1 / (gas.molecular_weights[gas.species_index('H')] * 2) + YO2Air / gas.molecular_weights[gas.species_index('O')]
LZC = 1 / KZ * (LZpCTable['H'] / (2 * gas.molecular_weights[gas.species_index('H')]) - LZpCTable['O'] / gas.molecular_weights[gas.species_index('O')])
for i in range(nPointsZ):
    gammaZYc[:, i] = LZC[:, i] / YkTable['H2O'][:, i][-1]
gammaZYc = gammaZYc * rhoTable

# Gamma_ZZ
gammaZZ      = np.zeros((nPointsC, nPointsZ))

if Soret:
    LambdakTable = {sp: m.calculateLambdaYkZ(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, ZTableExt1, species, sp, windowZ, order) for sp in species}
else:
    LambdakTable = {sp: m.calculateLambdaYkZ(DkTable, YkTable, MMWTable, ZTableExt1, species, sp, windowZ, order) for sp in species}

LZpZTable = {}
for k in listSpecies:
    isp = gas.element_index(k)
    S = sum(
        gas.n_atoms(gas.species_index(sp), isp) / gas.molecular_weights[gas.species_index(sp)] * LambdakTable[sp]
        for sp in species
    )
    LZpZTable[k] = gas.molecular_weights[gas.species_index(k)] * S

LZZ = 1 / KZ * (LZpZTable['H'] / (2 * gas.molecular_weights[gas.species_index('H')]) - LZpZTable['O'] / gas.molecular_weights[gas.species_index('O')])
for i in range(nPointsC):
    dcdz_YcConst  = -1 / YkTable['H2O'][-1, :] * (CTableExt[i, :] * dH2OdZ[-1, :])
    gammaZZ[i, :] = LZC[i, :] * dcdz_YcConst + LZZ[i, :] - DTable[i, :]
gammaZZ = (gammaZZ + DTable) * rhoTable


# #  Create G-HH
# gammaHH= np.zeros((nPointsC,nPointsZ))
# dTdH = m.SavitzkyGolay(YkTable['H2O'], hTable, 0, 40, order)

# LambdakTable = {}
# for k in species:
#     LambdakTable[k] = m.calculateLambdaYkh(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, hTable, species, k, 40, order)
    

# LHHTable = np.zeros((nPointsC,nPointsZ))
# for k in species:
#     LHHTable += LambdakTable[k]*hkTable[k]
    
# gammaHH = (kTable/rhoTable*dTdH + LHHTable)*rhoTable


# Extend Tables out of FL 
ZLean = np.linspace(0, zLean-1e-4, 50)
ZFlam = np.linspace(ZVec[0], ZVec[-1], 200)
ZRich = np.linspace(zRich+1e-4, 1, 50)

Z = np.concatenate((ZLean, ZFlam, ZRich))
ZTableExt = np.transpose(np.tile(Z[:, None], (1, nPointsC)))
CTableExt = np.tile(CVec[:, None], (1, np.size(Z)))

fsLVec = interp1d(np.concatenate(([0],ZVec,[1])), np.concatenate(([0],sLVec.flatten(),[0])), fill_value='extrapolate')
flFVec = interp1d(np.concatenate(([0],ZVec,[1])), np.concatenate(([0],lFVec.flatten(),[0])), fill_value='extrapolate')
ftauVec = interp1d(np.concatenate(([0],ZVec,[1])), np.concatenate(([0],tauVec.flatten(),[0])), fill_value='extrapolate')

sLVec = fsLVec(Z)
lFVec = flFVec(Z)
tauVec = ftauVec(Z)

TTableExt = np.zeros((nPointsC,len(Z)))
CpTableExt = np.zeros((nPointsC,len(Z)))
DTableExt = np.zeros((nPointsC,len(Z)))
dHTableExt = np.zeros((nPointsC,len(Z)))
kTableExt = np.zeros((nPointsC,len(Z)))
rhoTableExt = np.zeros((nPointsC,len(Z)))
MMWTableExt = np.zeros((nPointsC,len(Z)))
HRRTableExt = np.zeros((nPointsC,len(Z)))
sourceYcVTableExt = np.zeros((nPointsC,len(Z)))
gammaYcYcTableExt = np.zeros((nPointsC,len(Z)))
gammaYcZTableExt = np.zeros((nPointsC,len(Z)))
gammaZYcTableExt = np.zeros((nPointsC,len(Z)))
gammaZZTableExt = np.zeros((nPointsC,len(Z)))

YkTableExt = {}
XkTableExt = {}
DkTableExt = {}
DkSoretTableExt = {}
omegakTableExt = {}
for sp in species:
    YkTableExt[sp] = np.zeros((nPointsC,len(Z)))
    XkTableExt[sp] = np.zeros((nPointsC,len(Z)))
    DkTableExt[sp] = np.zeros((nPointsC,len(Z)))
    DkSoretTableExt[sp] = np.zeros((nPointsC,len(Z)))
    omegakTableExt[sp] = np.zeros((nPointsC,len(Z)))


# Interoplate over Z
for i in range(nPointsC):
    TTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], TTable[i,:], tin, tin)
    sourceYcVTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], sourceYcVTable[i,:], 0, 0)
    CpTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], CpTable[i,:], Air["CpAir"], Fuel["CpFuel"])
    DTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], DTable[i,:], Air["DAir"], Fuel["DFuel"])
    dHTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], dHTable[i,:], 0, 0)
    kTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], kTable[i,:], Air["kAir"], Fuel["kFuel"])
    rhoTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], rhoTable[i,:], Air["rhoAir"], Fuel["rhoFuel"])
    MMWTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], MMWTable[i,:], Air["WAir"], Fuel["WFuel"])
    HRRTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], HRRTable[i,:], 0, 0)
    gammaYcYcTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], gammaYcYc[i,:], 0, 0)
    gammaYcZTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], gammaYcZ[i,:], 0, 0)
    gammaZYcTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], gammaZYc[i,:], 0, 0)
    gammaZZTableExt[i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], gammaZZ[i,:], 0, 0)

    
    for sp in species:
        YkTableExt[sp][i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], YkTable[sp][i,:], Air[f"Y_{sp}Air"], Fuel[f"Y_{sp}Fuel"])
        XkTableExt[sp][i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], XkTable[sp][i,:], Air[f"X_{sp}Air"], Fuel[f"X_{sp}Fuel"])
        DkTableExt[sp][i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], DkTable[sp][i,:], Air[f"D_{sp}Air"], Fuel[f"D_{sp}Fuel"])
        omegakTableExt[sp][i,:] = m.interpExt2(ZTableExt[i,:], ZTableExt1[i,:], omegakTable[sp][i,:], 0, 0)


# 3D Plot of tables        
if plotTables == True:
    
    m.plot3DTable(CTableExt, ZTableExt, omegakTableExt['H2O'] ,"[kg/(m^3 \cdot s)]","\dot{\omega}_{\mathrm{H_2O}}", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, HRRTableExt ,"[W/m^3]","HRR", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, sourceYcVTableExt ,"[kg/(m^3 \cdot s)]","\dot{\omega}_{Yc,V}", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, TTableExt ,"[K]","T", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, DTableExt ,"[kg/(m \cdot s)]","D_{th}", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, CpTableExt ,"[J/(kg \cdot K)]","Cp", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, kTableExt ,"[W/(m \cdot K)]","\lambda", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, dHTableExt ,"[J/kg]","\Delta h_f^0", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, MMWTableExt ,"[kg/mol]","W_{mix}", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, YkTableExt['H2O'] ,"[-]","Y_{H2O}", "C", "Z")
    m.plot3DTable(CTableExt, ZTableExt, YkTableExt['OH'] ,"[-]","Y_{OH*}", "C", "Z")
    m.plot3DTable(ZTableExt, CTableExt, gammaYcYcTableExt ,"[m^2/s]","\\rho \cdot \Gamma_{Y_c,Y_c}", "Z", "C")
    m.plot3DTable(ZTableExt, CTableExt, gammaYcZTableExt ,"[m^2/s]","\\rho \cdot \Gamma_{Y_c,Z}", "Z", "C")
    m.plot3DTable(ZTableExt, CTableExt, gammaZYcTableExt ,"[m^2/s]","\\rho \cdot \Gamma'_{Z,Y_c}", "Z", "C") 
    m.plot3DTable(ZTableExt, CTableExt, gammaZZTableExt ,"[m^2/s]","\\rho \cdot \Gamma_{Z,Z}", "Z", "C") 
    
    
    # Contour Plot of tables
    m.plotContourTable(CTableExt, ZTableExt, HRRTableExt ,"[W/m^3]","HRR", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, omegakTableExt['H2O'] ,"[kg/(m^3 \cdot s)]","\dot{\omega}_{\mathrm{H_2O}}", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, sourceYcVTableExt ,"[kg/(m^3 \cdot s)]","\dot{\omega}_{Yc,V}", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, TTableExt ,"[K]","T", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, DTableExt ,"[kg/(m \cdot s)]","D_{th}", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, CpTableExt ,"[J/(kg \cdot K)]","Cp", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, dHTableExt ,"[J/kg]","\Delta h_f^0", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, kTableExt ,"[W/(m \cdot K)]","\lambda", "C", "Z")
    m.plotContourTable(CTableExt, ZTableExt, MMWTableExt,"[kg/mol]","W_{mix}", "C", "Z")
    m.plotContourTable(ZTableExt, CTableExt, gammaYcYcTableExt ,"[m^2/s]","\\rho \cdot \Gamma_{Y_c,Y_c}", "Z", "C")
    m.plotContourTable(ZTableExt, CTableExt, gammaYcZTableExt ,"[m^2/s]","\\rho \cdot \Gamma_{Y_c,Z}", "Z", "C")
    m.plotContourTable(ZTableExt, CTableExt, gammaZYcTableExt ,"[m^2/s]","\\rho \cdot \Gamma'_{Z,Y_c}", "Z", "C") 
    m.plotContourTable(ZTableExt, CTableExt, gammaZZTableExt ,"[m^2/s]","\\rho \cdot \Gamma_{Z,Z}", "Z", "C") 


# Save Tables
if save == True:
    
    os.makedirs(pathOut, exist_ok=True)

    gammaYcYc = pd.DataFrame(gammaYcYcTableExt, columns=Z)
    gammaYcZ = pd.DataFrame(gammaYcZTableExt, columns=Z)
    gammaZYc = pd.DataFrame(gammaZYcTableExt, columns=Z)
    gammaZZ = pd.DataFrame(gammaZZTableExt, columns=Z)
    CTable = pd.DataFrame(CTableExt, columns=Z)
    ZTable = pd.DataFrame(ZTableExt, columns=Z)
    sLVec = pd.DataFrame(sLVec, index=Z, columns=['sL'])
    tauVec = pd.DataFrame(tauVec, index=Z, columns=['tau'])
    lFVec = pd.DataFrame(lFVec, index=Z, columns=['lF'])
    TTable = pd.DataFrame(TTableExt, columns=Z)
    sourceYcVTable = pd.DataFrame(sourceYcVTableExt, columns=Z)
    DTable = pd.DataFrame(DTableExt, columns=Z)
    rhoTable = pd.DataFrame(rhoTableExt, columns=Z)
    dHTable = pd.DataFrame(dHTableExt, columns=Z)
    kTable = pd.DataFrame(kTableExt, columns=Z)
    CpTable = pd.DataFrame(CpTableExt, columns=Z)
    HRRTable = pd.DataFrame(HRRTableExt, columns=Z)
    MMWTable = pd.DataFrame(MMWTableExt, columns=Z)
    for i, sp in enumerate(species):
        omegakTable[sp] = pd.DataFrame(omegakTableExt[sp], columns=Z)
        YkTable[sp] = pd.DataFrame(YkTableExt[sp], columns=Z)
        XkTable[sp] = pd.DataFrame(XkTableExt[sp], columns=Z)
        DkTable[sp] = pd.DataFrame(DkTableExt[sp], columns=Z)
    
    gammaYcYc.to_csv(pathOut + "GYcYcTable.csv", index=False)
    gammaYcZ.to_csv(pathOut + "GYcZTable.csv", index=False)
    gammaZYc.to_csv(pathOut + "GZYcTable.csv", index=False)
    gammaZZ.to_csv(pathOut + "GZZTable.csv", index=False)
    CTable.to_csv(pathOut + "CTable.csv", index=False)
    sourceYcVTable.to_csv(pathOut + "sourceYcVTable.csv", index=False)
    ZTable.to_csv(pathOut + "ZTable.csv", index=False)
    sLVec.to_csv(pathOut + "sLVec.csv", index=False)
    lFVec.to_csv(pathOut + "lFVec.csv", index=False)
    tauVec.to_csv(pathOut + "tauVec.csv", index=False)
    TTable.to_csv(pathOut + "TTable.csv", index=False)
    dHTable.to_csv(pathOut + "dHTable.csv", index=False)
    DTable.to_csv(pathOut + "DTable.csv", index=False)
    rhoTable.to_csv(pathOut + "rhoTable.csv", index=False)
    kTable.to_csv(pathOut + "kTable.csv", index=False)
    CpTable.to_csv(pathOut + "CpTable.csv", index=False)
    HRRTable.to_csv(pathOut + "HRRTable.csv", index=False)
    MMWTable.to_csv(pathOut + "MMWTable.csv", index=False)
    for i, sp in enumerate(species):
        omegakTable[sp].to_csv(pathOut + f"omega_{sp}Table.csv", index=False)
        DkTable[sp].to_csv(pathOut + f"D_{sp}Table.csv", index=False)
        YkTable[sp].to_csv(pathOut + f"Y_{sp}Table.csv", index=False)
        XkTable[sp].to_csv(pathOut + f"X_{sp}Table.csv", index=False)
        
    print("Tables Succesfully saved in folder ---> " + nameMech + "Extended!")
    



    
    
