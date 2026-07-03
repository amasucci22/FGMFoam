# -*- coding: utf-8 -*-
"""
FGM TABLE GENERATION - Parallel version
Creates FGM tables X = X(C,Z) from freely-propagating premixed flat flames.
Each flamelet (phi) is solved independently and in parallel.
"""

import pandas as pd
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mpl
from multiprocessing import Pool, Value
import ctypes
import modules.functions as m

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern']

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
save        = True
plotTables  = False
Soret       = False

nameMech = 'SanDiegoNO'
mech     = nameMech + '.yaml'

if Soret:
    path     = "./" + nameMech + "Soret/"
    print("Soret activated!")
else:
    path     = "./" + nameMech + "/"
    print("Soret disactivated!")

CANTERA_DATA_DIR = "C:/Program Files/Cantera/data"

# Number of parallel workers (leave 1 core free)
N_WORKERS = max(1, os.cpu_count() - 1)

# ─────────────────────────────────────────────────────────────────────────────
# DISCRETIZATION
# ─────────────────────────────────────────────────────────────────────────────
nPointsZ  = 300
a         = 3
s         = np.linspace(0, 1, nPointsZ)
phiLean   = 0.28
phiRich   = 5.0
phiVector = phiLean + (phiRich - phiLean) * (np.exp(a * s) - 1) / (np.exp(a) - 1)

nPointsC  = 600
a         = -6  
x         = np.linspace(0, 1, nPointsC)
c         = (np.exp(a * x) - 1) / (np.exp(a) - 1)


# ─────────────────────────────────────────────────────────────────────────────
# WORKER FUNCTION  (runs in a separate process for each phi)
# ─────────────────────────────────────────────────────────────────────────────
# Shared counter initializer — called once per worker process at startup
def _init_counter(counter, total):
    global _counter, _total
    _counter = counter
    _total   = total


def runFlamelet(args):
    """
    Solve one premixed flamelet and return all tabulated data for that phi.
    Must be a top-level function so multiprocessing can pickle it.
    All imports are done inside the function to avoid pickling issues with
    Cantera objects.
    """
    i, phi, p, tin, mech, c, YScal, Soret, cantera_dir = args

    # Re-import inside worker (Cantera objects are not safely picklable)
    import cantera as ct
    import numpy as np
    import modules.functions as m

    ct.add_directory(cantera_dir)

    # ── Solve flamelet ────────────────────────────────────────────────────────
    f, gas = m.premixedFlame(p, tin, phi, mech, Soret=Soret)

    with _counter.get_lock():
        _counter.value += 1
        done = _counter.value
    pct = done / _total * 100
    print(f"  Flamelets completed: {done}/{_total}  ({pct:.1f}%)", flush=True)
    species = gas.species_names

    # Fresh gas object at this phi for mixture-fraction calculation
    gas2 = ct.Solution(mech)
    air  = "O2:0.21,N2:0.79"
    gas2.set_equivalence_ratio(phi=phi, fuel="H2:1", oxidizer=air)
    gas2.TPX = tin, p, gas2.X

    ZPhi = gas2.mixture_fraction(fuel="H2:1", oxidizer=air, element='Bilger')
    sL    = f.velocity[0]
    tau   = (f.T[-1] - f.T[0]) / f.T[0]
    lF    = (f.T[-1] - f.T[0]) / np.max(np.gradient(f.T, f.grid))

    # ── Bilger mixture fraction along the flame ───────────────────────────────
    ZBilg = m.BilgerMixtureFraction(gas2, f, mech)

    # ── Chemical enthalpy ─────────────────────────────────────────────────────
    hSpecies = (gas2.standard_enthalpies_RT * ct.gas_constant * gas2.T) / gas2.molecular_weights
    hChem    = np.dot(hSpecies, f.Y)

    # ── Scale progress variable Yc -> C ──────────────────────────────────────
    YcMax  = np.interp(ZBilg, YScal[:, 0], YScal[:, 1])
    cStar  = f.Y[gas2.species_index('H2O')] / YcMax
    cStar, idx = np.unique(cStar, return_index=True)
    cStar  = (cStar - cStar.min()) / (cStar.max() - cStar.min())

    # Helper: interpolate a 1-D array from cStar -> c
    def interp(arr):
        return m.interpExt(c, cStar, arr[idx])

    # ── Build data dictionary ─────────────────────────────────────────────────
    data = {
        "x"         : interp(f.grid),
        "T"         : interp(f.T),
        "U"         : interp(f.velocity),
        "rho"       : interp(f.density),
        "Z"         : interp(ZBilg),
        "Yc"        : interp(f.Y[gas2.species_index('H2O')]),
        "c"         : c,
        "Cp"        : interp(f.cp_mass),
        "k"         : interp(f.thermal_conductivity),
        "D"         : interp(f.thermal_conductivity / (f.density * f.cp_mass)),
        "h"         : interp(f.enthalpy_mass),
        "dH"        : interp(hChem),
        "HRR"       : interp(f.heat_release_rate),
        "MMW"       : interp(f.mean_molecular_weight) * 1e-3,
        # scalar flame properties
        "sL"        : sL,
        "tau"       : tau,
        "lF"        : lF,
        "ZPhi"     : ZPhi,
    }

    for sp in species:
        sp_idx = gas2.species_index(sp)
        data[f"D_{sp}"]      = interp(f.mix_diff_coeffs[sp_idx])
        data[f"DSoret_{sp}"] = interp(f.thermal_diff_coeffs[sp_idx])
        data[f"Y_{sp}"]      = interp(f.Y[sp_idx])
        data[f"X_{sp}"]      = interp(f.X[sp_idx])
        data[f"h_{sp}"]      = interp(f.partial_molar_enthalpies[sp_idx] / gas2.molecular_weights[sp_idx])
        data[f"omega_{sp}"]  = interp(f.net_production_rates[sp_idx] * gas2.molecular_weights[sp_idx])

    data["sourceYcV"] = interp(
        f.net_production_rates[gas2.species_index('H2O')]
        * gas2.molecular_weights[gas2.species_index('H2O')]
        * f.Y[gas2.species_index('H2O')]
    )

    return i, phi, data


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    ct.add_directory(CANTERA_DATA_DIR)

    # ── Load scaling file ─────────────────────────────────────────────────────
    YScal = pd.read_csv(
        path + "scaleYc" + nameMech + ".csv",
        header=None, skiprows=1
    ).to_numpy()

    # ── Air / Fuel boundary conditions ────────────────────────────────────────
    p   = ct.one_atm
    tin = 298.15
    
    # Create Air
    gasAir = ct.Solution(mech)
    species = gasAir.species_names
    air = "O2:0.21,N2:0.79"
    gasAir.set_equivalence_ratio(phi=0, fuel="H2:1", oxidizer=air)
    gasAir.TPX = tin, p, gasAir.X
    YO2Air = gasAir.Y[gasAir.species_index('O2')]

    # ── Build argument list for workers ──────────────────────────────────────
    args_list = [
        (i, phi, p, tin, mech, c, YScal, Soret, CANTERA_DATA_DIR)
        for i, phi in enumerate(phiVector)
    ]

    # ── Run flamelets in parallel ─────────────────────────────────────────────
    print(f"\nRunning {nPointsZ} flamelets on {N_WORKERS} workers...\n")

    counter = Value(ctypes.c_int, 0)
    with Pool(processes=N_WORKERS,
              initializer=_init_counter,
              initargs=(counter, nPointsZ)) as pool:
        rawResults = pool.map(runFlamelet, args_list)

    # Sort by original phi index to preserve order
    rawResults.sort(key=lambda x: x[0])

    # Rebuild results dict and scalar vectors
    results  = {}
    ZVec     = np.zeros(nPointsZ)
    sLVec    = np.zeros(nPointsZ)
    tauVec   = np.zeros(nPointsZ)
    lFVec    = np.zeros(nPointsZ)

    for i, phi, data in rawResults:
        results[phi]  = data
        ZVec[i]       = data["ZPhi"]
        sLVec[i]      = data["sL"]
        tauVec[i]     = data["tau"]
        lFVec[i]      = data["lF"]

    print("\nAll flamelets done. Building tables...")

    # ── Build 2-D tables (nPointsC x nPointsZ) ───────────────────────────────
    def makeTable(key):
        return m.createTables(phiVector, key, results, nPointsC)

    CTable         = makeTable("c");        CTable[CTable < 1e-6] = 0
    ZTable         = makeTable("Z")
    TTable         = makeTable("T")
    DTable         = makeTable("D")
    rhoTable       = makeTable("rho")
    dHTable        = makeTable("dH")
    hTable         = makeTable("h")
    kTable         = makeTable("k")
    CpTable        = makeTable("Cp")
    HRRTable       = makeTable("HRR")
    MMWTable       = makeTable("MMW")
    sourceYcVTable = makeTable("sourceYcV")

    omegakTable    = {sp: makeTable(f"omega_{sp}") for sp in species}
    YkTable        = {sp: makeTable(f"Y_{sp}")     for sp in species}
    hkTable        = {sp: makeTable(f"h_{sp}")     for sp in species}
    XkTable        = {sp: makeTable(f"X_{sp}")     for sp in species}
    DkTable        = {sp: makeTable(f"D_{sp}")     for sp in species}
    DkSoretTable   = {sp: makeTable(f"DSoret_{sp}") for sp in species}

    # ── Convert to DataFrames and optionally save ─────────────────────────────
    def to_df(arr):
        return pd.DataFrame(arr, columns=phiVector)

    CTable_df         = to_df(CTable)
    ZTable_df         = to_df(ZTable)
    TTable_df         = to_df(TTable)
    DTable_df         = to_df(DTable)
    rhoTable_df       = to_df(rhoTable)
    dHTable_df        = to_df(dHTable)
    hTable_df         = to_df(hTable)
    kTable_df         = to_df(kTable)
    CpTable_df        = to_df(CpTable)
    HRRTable_df       = to_df(HRRTable)
    MMWTable_df       = to_df(MMWTable)
    sourceYcVTable_df = to_df(sourceYcVTable)

    sLVec_df  = pd.DataFrame(sLVec,  index=phiVector, columns=['sL'])
    tauVec_df = pd.DataFrame(tauVec, index=phiVector, columns=['tau'])
    lFVec_df  = pd.DataFrame(lFVec,  index=phiVector, columns=['lF'])

    omegakTable_df  = {sp: to_df(omegakTable[sp])  for sp in species}
    YkTable_df      = {sp: to_df(YkTable[sp])      for sp in species}
    hkTable_df      = {sp: to_df(hkTable[sp])      for sp in species}
    XkTable_df      = {sp: to_df(XkTable[sp])      for sp in species}
    DkTable_df      = {sp: to_df(DkTable[sp])      for sp in species}
    DkSoretTable_df = {sp: to_df(DkSoretTable[sp]) for sp in species}

    if save:
        CTable_df.to_csv(path + "CTable.csv",             index=False)
        sLVec_df.to_csv( path + "sLVec.csv",              index=False)
        lFVec_df.to_csv( path + "lFVec.csv",              index=False)
        tauVec_df.to_csv(path + "tauVec.csv",             index=False)
        ZTable_df.to_csv(path + "ZTable.csv",             index=False)
        TTable_df.to_csv(path + "TTable.csv",             index=False)
        DTable_df.to_csv(path + "DTable.csv",             index=False)
        rhoTable_df.to_csv(path + "rhoTable.csv",         index=False)
        dHTable_df.to_csv(path + "dHTable.csv",           index=False)
        hTable_df.to_csv(path + "hTable.csv",             index=False)
        kTable_df.to_csv(path + "kTable.csv",   index=False)
        CpTable_df.to_csv(path + "CpTable.csv",           index=False)
        HRRTable_df.to_csv(path + "HRRTable.csv",         index=False)
        MMWTable_df.to_csv(path + "MMWTable.csv",         index=False)
        sourceYcVTable_df.to_csv(path + "sourceYcVTable.csv", index=False)
        for sp in species:
            omegakTable_df[sp].to_csv( path + f"omega_{sp}Table.csv",  index=False)
            DkTable_df[sp].to_csv(     path + f"D_{sp}Table.csv",      index=False)
            DkSoretTable_df[sp].to_csv(path + f"DSoret_{sp}Table.csv", index=False)
            YkTable_df[sp].to_csv(     path + f"Y_{sp}Table.csv",      index=False)
            hkTable_df[sp].to_csv(     path + f"h_{sp}Table.csv",      index=False)
            XkTable_df[sp].to_csv(     path + f"X_{sp}Table.csv",      index=False)
        print("Tables successfully saved in folder ---> " + nameMech + "!")

    # ── Back to numpy for post-processing ────────────────────────────────────
    CTable         = CTable_df.to_numpy()
    ZTable         = ZTable_df.to_numpy()
    TTable         = TTable_df.to_numpy()
    DTable         = DTable_df.to_numpy()
    hTable         = hTable_df.to_numpy()
    tauVec         = tauVec_df.to_numpy()
    sLVec          = sLVec_df.to_numpy()
    lFVec          = lFVec_df.to_numpy()
    rhoTable       = rhoTable_df.to_numpy()
    dHTable        = dHTable_df.to_numpy()
    kTable         = kTable_df.to_numpy()
    CpTable        = CpTable_df.to_numpy()
    HRRTable       = HRRTable_df.to_numpy()
    MMWTable       = MMWTable_df.to_numpy()
    sourceYcVTable = sourceYcVTable_df.to_numpy()
    for sp in species:
        omegakTable[sp]  = omegakTable_df[sp].to_numpy()
        YkTable[sp]      = YkTable_df[sp].to_numpy()
        hkTable[sp]      = hkTable_df[sp].to_numpy()
        XkTable[sp]      = XkTable_df[sp].to_numpy()
        DkTable[sp]      = DkTable_df[sp].to_numpy()
        DkSoretTable[sp] = DkSoretTable_df[sp].to_numpy()

    # ── Create Z/C grids ────────────────────────────────────────────────────
    ZTableExt = np.zeros((nPointsC, nPointsZ))
    CTableExt = CTable
    for i in range(nPointsZ):
        ZTableExt[:, i] = ZTable[0, i]

    # ── Interpolate all tables over Z ─────────────────────────────────────────
    for i in range(nPointsC):
        ZiSRC = ZTable[i, :]
        ZiDST = ZTableExt[i, :]
        TTable[i, :]         = m.interpExt(ZiDST, ZiSRC, TTable[i, :])
        CpTable[i, :]        = m.interpExt(ZiDST, ZiSRC, CpTable[i, :])
        DTable[i, :]         = m.interpExt(ZiDST, ZiSRC, DTable[i, :])
        kTable[i, :]         = m.interpExt(ZiDST, ZiSRC, kTable[i, :])
        rhoTable[i, :]       = m.interpExt(ZiDST, ZiSRC, rhoTable[i, :])
        dHTable[i, :]        = m.interpExt(ZiDST, ZiSRC, dHTable[i, :])
        hTable[i, :]         = m.interpExt(ZiDST, ZiSRC, hTable[i, :])
        MMWTable[i, :]       = m.interpExt(ZiDST, ZiSRC, MMWTable[i, :])
        HRRTable[i, :]       = m.interpExt(ZiDST, ZiSRC, HRRTable[i, :])
        sourceYcVTable[i, :] = m.interpExt(ZiDST, ZiSRC, sourceYcVTable[i, :])
        for sp in species:
            YkTable[sp][i, :]      = m.interpExt(ZiDST, ZiSRC, YkTable[sp][i, :])
            hkTable[sp][i, :]      = m.interpExt(ZiDST, ZiSRC, hkTable[sp][i, :])
            XkTable[sp][i, :]      = m.interpExt(ZiDST, ZiSRC, XkTable[sp][i, :])
            DkTable[sp][i, :]      = m.interpExt(ZiDST, ZiSRC, DkTable[sp][i, :])
            DkSoretTable[sp][i, :] = m.interpExt(ZiDST, ZiSRC, DkSoretTable[sp][i, :])
            omegakTable[sp][i, :]  = m.interpExt(ZiDST, ZiSRC, omegakTable[sp][i, :])

    # ── Gamma tensors ─────────────────────────────────────────────────────────
    windowC = 40
    windowZ = 40
    order   = 1

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
    dH2OdZ    = m.SavitzkyGolay(YkTable['H2O'], ZTableExt, 1, windowC, order)
    
    if Soret:
        LYcZ  = m.calculateLambdaYkZ(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, ZTableExt, species, 'H2O', windowC, order)
    else:
        LYcZ  = m.calculateLambdaYkZ(DkTable, YkTable, MMWTable, ZTableExt, species, 'H2O', windowC, order)

    for i in range(nPointsC):
        dcdz_YcConst  = -1 / YkTable['H2O'][-1, :] * (CTableExt[i, :] * dH2OdZ[-1, :])
        gammaYcZ[i,:] = LYcC[i, :] * dcdz_YcConst + LYcZ[i, :]
    gammaYcZ = gammaYcZ * rhoTable

    # Gamma_ZYc
    gammaZYc = np.zeros((nPointsC, nPointsZ))
    
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
        LambdakTable = {sp: m.calculateLambdaYkZ(DkTable, DkSoretTable, TTable, rhoTable, YkTable, MMWTable, ZTableExt, species, sp, windowZ, order) for sp in species}
    else:
        LambdakTable = {sp: m.calculateLambdaYkZ(DkTable, YkTable, MMWTable, ZTableExt, species, sp, windowZ, order) for sp in species}
  
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

    # ── Clip source terms ─────────────────────────────────────────────────────
    omegakTable['H2O'][omegakTable['H2O'] < 0]    = 0
    omegakTable['H2O'][omegakTable['H2O'] < 1e-10] = 0
    omegakTable['H2'][omegakTable['H2'] > 0]       = 0
    omegakTable['H2'][omegakTable['H2'] > 1e-10]   = 0

    ZVec = ZTableExt[0, :]
    CVec = CTableExt[:, 0]

    # ── 3-D / Contour plots ───────────────────────────────────────────────────
    if plotTables:
        m.plot3DTable(CTableExt, ZTableExt, omegakTable['H2O'], "[kg/(m^3 \\cdot s)]", "\\dot{\\omega}_{\\mathrm{H_2O}}", "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, HRRTable,           "[W/m^3]",              "HRR",                              "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, TTable,             "[K]",                  "T",                                "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, DTable,             "[kg/(m \\cdot s)]",    "D_{th}",                           "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, DkSoretTable['H2O'],"[m^2/s]",              "DSoret_{H_2O}",                    "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, DkTable['H2O'],     "[m^2/s]",              "D_{H_2O}",                         "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, CpTable,            "[J/(kg \\cdot K)]",    "Cp",                               "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, kTable,        "[W/(m \\cdot K)]",     "\\lambda",                         "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, dHTable,            "[J/kg]",               "\\Delta h_f^0",                    "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, MMWTable,           "[kg/mol]",             "W_{mix}",                          "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, YkTable['H2O'],     "[-]",                  "Y_{H2O}",                          "C", "Z")
        m.plot3DTable(CTableExt, ZTableExt, YkTable['OH'],      "[-]",                  "Y_{OH*}",                          "C", "Z")
        m.plot3DTable(ZTableExt, CTableExt, gammaYcYc,          "[m^2/s]", "\\rho \\cdot \\Gamma'_{Y_c,Yc}",               "Z", "C")
        m.plot3DTable(ZTableExt, CTableExt, gammaYcZ,           "[m^2/s]", "\\rho \\cdot \\Gamma'_{Y_c,Z}",                "Z", "C")
        m.plot3DTable(ZTableExt, CTableExt, gammaZYc,           "[m^2/s]", "\\rho \\cdot \\Gamma'_{Z,Y_c}",                "Z", "C")
        m.plot3DTable(ZTableExt, CTableExt, gammaZZ,            "[m^2/s]", "\\rho \\cdot \\Gamma_{Z,Z}",                   "Z", "C")

        m.plotContourTable(CTableExt, ZTableExt, HRRTable,            "[W/m^3]",           "HRR",                            "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, omegakTable['H2O'],  "[kg/(m^3 \\cdot s)]","\\dot{\\omega}_{\\mathrm{H_2O}}","C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, TTable,              "[K]",               "T",                              "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, DTable,              "[kg/(m \\cdot s)]", "D_{th}",                         "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, DkSoretTable['H2O'], "[m^2/s]",           "DSoret_{H_2O}",                  "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, DkTable['H2O'],      "[m^2/s]",           "D_{H_2O}",                       "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, CpTable,             "[J/(kg \\cdot K)]", "Cp",                             "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, dHTable,             "[J/kg]",            "\\Delta h_f^0",                  "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, kTable,         "[W/(m \\cdot K)]",  "\\lambda",                       "C", "Z")
        m.plotContourTable(CTableExt, ZTableExt, MMWTable,            "[kg/mol]",          "W_{mix}",                        "C", "Z")
        m.plotContourTable(ZTableExt, CTableExt, gammaYcYc,           "[m^2/s]", "\\rho \\cdot \\Gamma_{Y_c,Y_c}",          "Z", "C")
        m.plotContourTable(ZTableExt, CTableExt, gammaYcZ,            "[m^2/s]", "\\rho \\cdot \\Gamma_{Y_c,Z}",            "Z", "C")
        m.plotContourTable(ZTableExt, CTableExt, gammaZYc,            "[m^2/s]", "\\rho \\cdot \\Gamma'_{Z,Y_c}",           "Z", "C")
        m.plotContourTable(ZTableExt, CTableExt, gammaZZ,             "[m^2/s]", "\\rho \\cdot \\Gamma_{Z,Z}",              "Z", "C")

    print("\nDone.")
