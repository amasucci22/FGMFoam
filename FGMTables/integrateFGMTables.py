# -*- coding: utf-8 -*-
"""
INTEGRATE FGM TABLES X = X(C,Z) with Beta-PDF  —  Parallel version

"""

import cantera as ct
import pandas as pd
import numpy as np
import matplotlib as mpl
from scipy.stats import beta
import os
import ctypes
from multiprocessing import Pool, Value
import modules.functions as m


# ─────────────────────────────────────────────────────────────────────────────
# PROGRESS BAR  (terminal, no external dependencies)
# ─────────────────────────────────────────────────────────────────────────────
def _print_bar(done, total, label, bar_width=40):
    frac    = done / total
    filled  = int(bar_width * frac)
    bar     = '█' * filled + '░' * (bar_width - filled)
    pct     = frac * 100
    end     = '\n' if done == total else '\r'

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
save        = True
plotTables  = False
exportTable = True
NOEqn       = True
Soret       = False

nameMech = 'SanDiegoNO'
mech     = nameMech + '.yaml'

N_WORKERS = max(1, os.cpu_count() - 1)

if Soret:
    path  = "./" + nameMech + "SoretExtended/"
    YScal = pd.read_csv(nameMech + "Soret/scaleYc" + nameMech + ".csv", header=None, skiprows=1).to_numpy()
    ZTable = pd.read_csv(nameMech + "Soret/ZTable.csv", header=None, skiprows=1).to_numpy()
    print("Soret activated!")
else:
    path  = "./" + nameMech + "Extended/"
    YScal = pd.read_csv(nameMech + "/scaleYc" + nameMech + ".csv", header=None, skiprows=1).to_numpy()
    ZTable = pd.read_csv(nameMech + "/ZTable.csv", header=None, skiprows=1).to_numpy()
    print("Soret disactivated!")

# ─────────────────────────────────────────────────────────────────────────────
# DISCRETIZATION
# ─────────────────────────────────────────────────────────────────────────────
nCTilde = 150
nZTilde = 120
nCVar   = 11
nZVar   = 3

nPointsZLean = 15
nPointsZRich = 15

ZVec  = ZTable[0, :]
zLean = ZVec[0]
zRich = ZVec[-1]

a = 3
s = np.linspace(0, 1, nZTilde)
ZVec = zLean + (zRich - zLean) * (np.exp(a * s) - 1) / (np.exp(a) - 1)

# Z lean / rich extensions
ZLean  = np.linspace(0, zLean - 1e-4, nPointsZLean)
ZRich  = np.linspace(zRich + 1e-4, 1, nPointsZRich)
ZTilde01 = np.concatenate((ZLean, ZVec, ZRich))
nZTilde  = len(ZTilde01)

a = -6
x = np.linspace(0, 1, nCTilde)
CTilde01 = (np.exp(a * x) - 1) / (np.exp(a) - 1)

p   = ct.one_atm
tin = 298.15

# ─────────────────────────────────────────────────────────────────────────────
# Air / Fuel boundary properties
# ─────────────────────────────────────────────────────────────────────────────
gas    = ct.Solution(mech)
gasAir = ct.Solution(mech)
species = gasAir.species_names
air     = "O2:0.21,N2:0.79"
gasAir.set_equivalence_ratio(phi=0, fuel="H2:1", oxidizer=air)
gasAir.TPX = tin, p, gasAir.X
YO2Air = gasAir.Y[gasAir.species_index('O2')]

Air = {
    "CpAir"     : gasAir.cp_mass,
    "WAir"      : gasAir.mean_molecular_weight * 1e-3,
    "DAir"      : gasAir.thermal_conductivity / (gasAir.density * gasAir.cp_mass),
    "rhoAir"    : gasAir.density,
    "kAir"      : gasAir.thermal_conductivity,
    "HRRAir"    : 0,
}
for sp in species:
    Air[f"D_{sp}Air"] = gasAir.mix_diff_coeffs_mass[gasAir.species_index(sp)]
    Air[f"Y_{sp}Air"] = gasAir.Y[gasAir.species_index(sp)]
    Air[f"X_{sp}Air"] = gasAir.X[gasAir.species_index(sp)]

gasFuel = ct.Solution(mech)
gasFuel.TPX = tin, p, "H2:1"
Fuel = {
    "CpFuel"     : gasFuel.cp_mass,
    "WFuel"      : gasFuel.mean_molecular_weight * 1e-3,
    "DFuel"      : gasFuel.thermal_conductivity / (gasFuel.density * gasFuel.cp_mass),
    "rhoFuel"    : gasFuel.density,
    "kFuel"      : gasFuel.thermal_conductivity,
    "HRRFuel"    : 0,
}
for sp in species:
    Fuel[f"D_{sp}Fuel"] = gasFuel.mix_diff_coeffs_mass[gasFuel.species_index(sp)]
    Fuel[f"Y_{sp}Fuel"] = gasFuel.Y[gasFuel.species_index(sp)]
    Fuel[f"X_{sp}Fuel"] = gasFuel.X[gasFuel.species_index(sp)]

# ─────────────────────────────────────────────────────────────────────────────
# LOAD FGM TABLES 
# ─────────────────────────────────────────────────────────────────────────────
DkTable     = {}
YkTable     = {}
omegakTable = {}
XkTable     = {}

CTable         = pd.read_csv(path + "CTable.csv",         header=None, skiprows=1).to_numpy()
ZTable         = pd.read_csv(path + "ZTable.csv",         header=None, skiprows=1).to_numpy()
sLVec          = pd.read_csv(path + "sLVec.csv",          header=None, skiprows=1).to_numpy()
tauVec         = pd.read_csv(path + "tauVec.csv",         header=None, skiprows=1).to_numpy()
lFVec          = pd.read_csv(path + "lFVec.csv",          header=None, skiprows=1).to_numpy()
TTable         = pd.read_csv(path + "TTable.csv",         header=None, skiprows=1).to_numpy()
DTable         = pd.read_csv(path + "DTable.csv",         header=None, skiprows=1).to_numpy()
rhoTable       = pd.read_csv(path + "rhoTable.csv",       header=None, skiprows=1).to_numpy()
dHTable        = pd.read_csv(path + "dHTable.csv",        header=None, skiprows=1).to_numpy()
kTable         = pd.read_csv(path + "kTable.csv",         header=None, skiprows=1).to_numpy()
CpTable        = pd.read_csv(path + "CpTable.csv",        header=None, skiprows=1).to_numpy()
MMWTable       = pd.read_csv(path + "MMWTable.csv",       header=None, skiprows=1).to_numpy()
HRRTable       = pd.read_csv(path + "HRRTable.csv",       header=None, skiprows=1).to_numpy()
GYcYcTable     = pd.read_csv(path + "GYcYcTable.csv",     header=None, skiprows=1).to_numpy()
GYcZTable      = pd.read_csv(path + "GYcZTable.csv",      header=None, skiprows=1).to_numpy()
GZYcTable      = pd.read_csv(path + "GZYcTable.csv",      header=None, skiprows=1).to_numpy()
GZZTable       = pd.read_csv(path + "GZZTable.csv",       header=None, skiprows=1).to_numpy()
sourceYcVTable = pd.read_csv(path + "sourceYcVTable.csv", header=None, skiprows=1).to_numpy()

for sp in species:
    omegakTable[sp] = pd.read_csv(path + f"omega_{sp}Table.csv", header=None, skiprows=1).to_numpy()
    DkTable[sp]     = pd.read_csv(path + f"D_{sp}Table.csv",     header=None, skiprows=1).to_numpy()
    YkTable[sp]     = pd.read_csv(path + f"Y_{sp}Table.csv",     header=None, skiprows=1).to_numpy()

nPointsZ = len(ZTable[0,:])
nPointsC = len(CTable[:,0])


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Interpolate tables onto CTilde / ZTilde grids
# ─────────────────────────────────────────────────────────────────────────────
# -- interpolate over C axis --------------------------------------------------
def _interp_col_C(arr):
    out = np.zeros((nCTilde, nPointsZ))
    for i in range(nPointsZ):
        out[:, i] = m.interpExt(CTilde01, CTable[:, 0], arr[:, i])
    return out

TCTilde        = _interp_col_C(TTable)
WCTilde        = _interp_col_C(MMWTable)
DCTilde        = _interp_col_C(DTable)
kCTilde        = _interp_col_C(kTable)
CpCTilde       = _interp_col_C(CpTable)
dHCTilde       = _interp_col_C(dHTable)
rhoCTilde      = _interp_col_C(rhoTable)
HRRCTilde      = _interp_col_C(HRRTable)
GYcYcCTilde    = _interp_col_C(GYcYcTable)
GYcZCTilde     = _interp_col_C(GYcZTable)
GZYcCTilde     = _interp_col_C(GZYcTable)
GZZCTilde      = _interp_col_C(GZZTable)
sourceYcCTilde = _interp_col_C(omegakTable['H2O'])
sourceH2CTilde = _interp_col_C(omegakTable['H2'])
sourceYcVCTilde= _interp_col_C(sourceYcVTable)
if NOEqn:
    sourceNOCTilde = _interp_col_C(omegakTable['NO'])
YkCTilde = {sp: _interp_col_C(YkTable[sp]) for sp in species}
DkCTilde = {sp: _interp_col_C(DkTable[sp]) for sp in species}

CCTilde = np.tile(CTilde01[:, None], (1, nPointsZ))


# -- interpolate over Z axis --------------------------------------------------
def _interp_row_Z(arr_CTilde):
    out = np.zeros((nCTilde, nZTilde))
    for i in range(nCTilde):
        out[i, :] = m.interpExt(ZTilde01, ZTable[0, :], arr_CTilde[i, :])
    return out

WZTilde        = _interp_row_Z(WCTilde)
DZTilde        = _interp_row_Z(DCTilde)
dHZTilde       = _interp_row_Z(dHCTilde)
kZTilde        = _interp_row_Z(kCTilde)
CpZTilde       = _interp_row_Z(CpCTilde)
rhoZTilde      = _interp_row_Z(rhoCTilde)
HRRZTilde      = _interp_row_Z(HRRCTilde)
TZTilde        = _interp_row_Z(TCTilde)
GYcYcZTilde    = _interp_row_Z(GYcYcCTilde)
GYcZZTilde     = _interp_row_Z(GYcZCTilde)
GZYcZTilde     = _interp_row_Z(GZYcCTilde)
GZZZTilde      = _interp_row_Z(GZZCTilde)
sourceYcZTilde = _interp_row_Z(sourceYcCTilde)
sourceH2ZTilde = _interp_row_Z(sourceH2CTilde)
sourceYcVZTilde= _interp_row_Z(sourceYcVCTilde)
if NOEqn:
    sourceNOZTilde = _interp_row_Z(sourceNOCTilde)
YkZTilde = {sp: _interp_row_Z(YkCTilde[sp]) for sp in species}
DkZTilde = {sp: _interp_row_Z(DkCTilde[sp]) for sp in species}

CZTilde = np.tile(CTilde01[:, None], (1, nZTilde))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Beta parameters for C integration
# ─────────────────────────────────────────────────────────────────────────────
iexp     = np.arange(1, nCVar + 1, 1)
CVarMax01 = 0.25 * ((iexp - 1) / (nCVar - 1)) ** 2.5

[CVAR01, CTILDE01] = np.meshgrid(CVarMax01, CTilde01)
CTILDE = CTILDE01[1:-1, 1:-1]
CVAR   = CVAR01[1:-1, 1:-1]
CVAR   = np.minimum(CTILDE * (1 - CTILDE), CVAR)
idx    = CVAR == CTILDE * (1 - CTILDE)
CVAR[idx] = CVAR[idx] - 1e-8

a_C = CTILDE * (CTILDE * (1 - CTILDE) / CVAR - 1)
b_C = a_C * (1 / CTILDE - 1)

nPOINTSC = 1000
C        = np.linspace(0, 1, nPOINTSC)
Z        = ZTilde01
[ZZ, CC] = np.meshgrid(Z, C)

# -- interpolate tables onto fine C / Z grids ---------------------------------
def _interp_fine_Z(arr_ZTilde):
    out = np.zeros((nPOINTSC, nZTilde))
    for i in range(nZTilde):
        out[:, i] = m.interpExt(CC[:, 0], CZTilde[:, 0], arr_ZTilde[:, i])
    return out

W_fine        = _interp_fine_Z(WZTilde)
T_fine        = _interp_fine_Z(TZTilde)
dH_fine       = _interp_fine_Z(dHZTilde)
D_fine        = _interp_fine_Z(DZTilde)
k_fine        = _interp_fine_Z(kZTilde)
Cp_fine       = _interp_fine_Z(CpZTilde)
rho_fine      = _interp_fine_Z(rhoZTilde)
HRR_fine      = _interp_fine_Z(HRRZTilde)
GYcYc_fine    = _interp_fine_Z(GYcYcZTilde)
GYcZ_fine     = _interp_fine_Z(GYcZZTilde)
GZYc_fine     = _interp_fine_Z(GZYcZTilde)
GZZ_fine      = _interp_fine_Z(GZZZTilde)
sourceYc_fine = _interp_fine_Z(sourceYcZTilde)
sourceH2_fine = _interp_fine_Z(sourceH2ZTilde)
sourceYcV_fine= _interp_fine_Z(sourceYcVZTilde)
if NOEqn:
    sourceNO_fine = _interp_fine_Z(sourceNOZTilde)
Yk_fine = {sp: _interp_fine_Z(YkZTilde[sp]) for sp in species}
Dk_fine = {sp: _interp_fine_Z(DkZTilde[sp]) for sp in species}

# -- Savitzky-Golay derivatives over C ----------------------------------------
dWdC        = m.SavitzkyGolay(W_fine,        CC, 0, 2, 1)
dTdC        = m.SavitzkyGolay(T_fine,        CC, 0, 2, 1)
ddHdC       = m.SavitzkyGolay(dH_fine,       CC, 0, 2, 1)
dDdC        = m.SavitzkyGolay(D_fine,        CC, 0, 2, 1)
dkdC        = m.SavitzkyGolay(k_fine,        CC, 0, 2, 1)
dCpdC       = m.SavitzkyGolay(Cp_fine,       CC, 0, 2, 1)
drhodC      = m.SavitzkyGolay(rho_fine,      CC, 0, 2, 1)
dHRRdC      = m.SavitzkyGolay(HRR_fine,      CC, 0, 2, 1)
dGYcYcdC    = m.SavitzkyGolay(GYcYc_fine,    CC, 0, 2, 1)
dGYcZdC     = m.SavitzkyGolay(GYcZ_fine,     CC, 0, 2, 1)
dGZYcdC     = m.SavitzkyGolay(GZYc_fine,     CC, 0, 2, 1)
dGZZdC      = m.SavitzkyGolay(GZZ_fine,      CC, 0, 2, 1)
dsourceYcdC = m.SavitzkyGolay(sourceYc_fine, CC, 0, 2, 1)
dsourceYcVdC= m.SavitzkyGolay(sourceYcV_fine,CC, 0, 2, 1)
dsourceH2dC = m.SavitzkyGolay(sourceH2_fine, CC, 0, 2, 1)
if NOEqn:
    dsourceNOdC = m.SavitzkyGolay(sourceNO_fine, CC, 0, 2, 1)
dYkdC = {sp: m.SavitzkyGolay(Yk_fine[sp], CC, 0, 2, 1) for sp in species}
dDkdC = {sp: m.SavitzkyGolay(Dk_fine[sp], CC, 0, 2, 1) for sp in species}


# ─────────────────────────────────────────────────────────────────────────────
# PARALLEL BETA-PDF INTEGRATION OVER C
# ─────────────────────────────────────────────────────────────────────────────


def _init_intC(shared_arrays, sp_list, noe, counter, total):
    global _arrC, _species_list, _NOEqn, _counter, _total
    _arrC         = shared_arrays
    _species_list = sp_list
    _NOEqn        = noe
    _counter      = counter
    _total        = total


def _worker_intC(args):
    """Compute the Beta-PDF integral over C for one (i, j) point."""
    i, j, a1, b1 = args

    import numpy as np
    from scipy.stats import beta as beta_dist

    with _counter.get_lock():
        _counter.value += 1
        done = _counter.value
    pct    = done / _total * 100
    filled = int(40 * done / _total)
    bar    = '█' * filled + '░' * (40 - filled)
    end    = '\n' if done == _total else '\r'
    print(f"  Integration over C   {done:>5}/{_total}  {pct:5.1f}%  [{bar}]", end=end, flush=True)

    C_loc = _arrC['C']
    CDF   = beta_dist.cdf(C_loc, a1, b1)

    result = {
        'i': i, 'j': j,
        'W':        CDF[-1] * _arrC['W_fine'][-1,:]         - np.trapezoid(CDF[:, None] * _arrC['dWdC'],        C_loc, axis=0),
        'dH':       CDF[-1] * _arrC['dH_fine'][-1,:]        - np.trapezoid(CDF[:, None] * _arrC['ddHdC'],       C_loc, axis=0),
        'D':        CDF[-1] * _arrC['D_fine'][-1,:]         - np.trapezoid(CDF[:, None] * _arrC['dDdC'],        C_loc, axis=0),
        'rho':      CDF[-1] * _arrC['rho_fine'][-1,:]       - np.trapezoid(CDF[:, None] * _arrC['drhodC'],      C_loc, axis=0),
        'k':       CDF[-1] * _arrC['k_fine'][-1,:]          - np.trapezoid(CDF[:, None] * _arrC['dkdC'],   C_loc, axis=0),
        'Cp':       CDF[-1] * _arrC['Cp_fine'][-1,:]        - np.trapezoid(CDF[:, None] * _arrC['dCpdC'],       C_loc, axis=0),
        'HRR':      CDF[-1] * 0                             - np.trapezoid(CDF[:, None] * _arrC['dHRRdC'],      C_loc, axis=0),
        'T':        CDF[-1] * _arrC['T_fine'][-1,:]         - np.trapezoid(CDF[:, None] * _arrC['dTdC'],        C_loc, axis=0),
        'GYcYc':    CDF[-1] * _arrC['GYcYc_fine'][-1,:]     - np.trapezoid(CDF[:, None] * _arrC['dGYcYcdC'],    C_loc, axis=0),
        'GYcZ':     CDF[-1] * _arrC['GYcZ_fine'][-1,:]      - np.trapezoid(CDF[:, None] * _arrC['dGYcZdC'],     C_loc, axis=0),
        'GZYc':     CDF[-1] * _arrC['GZYc_fine'][-1,:]      - np.trapezoid(CDF[:, None] * _arrC['dGZYcdC'],     C_loc, axis=0),
        'GZZ':      CDF[-1] * _arrC['GZZ_fine'][-1,:]       - np.trapezoid(CDF[:, None] * _arrC['dGZZdC'],      C_loc, axis=0),
        'sourceH2': CDF[-1] * 0                             - np.trapezoid(CDF[:, None] * _arrC['dsourceH2dC'], C_loc, axis=0),
        'sourceYc': CDF[-1] * 0                             - np.trapezoid(CDF[:, None] * _arrC['dsourceYcdC'], C_loc, axis=0),
        'sourceYcV':CDF[-1] * 0                             - np.trapezoid(CDF[:, None] * _arrC['dsourceYcVdC'],C_loc, axis=0),
    }
    if _NOEqn:
        result['sourceNO'] = (CDF[-1] * _arrC['sourceNO_fine'][-1,:]
                              - np.trapezoid(CDF[:, None] * _arrC['dsourceNOdC'], C_loc, axis=0))
    result['Yk'] = {sp: CDF[-1] * _arrC['Yk_fine'][sp][-1, :]
                        - np.trapezoid(CDF[:, None] * _arrC['dYkdC'][sp], C_loc, axis=0)
                    for sp in _species_list}
    result['Dk'] = {sp: CDF[-1] * _arrC['Dk_fine'][sp][-1, :]
                        - np.trapezoid(CDF[:, None] * _arrC['dDkdC'][sp], C_loc, axis=0)
                    for sp in _species_list}
    return result


# ─────────────────────────────────────────────────────────────────────────────
# PARALLEL BETA-PDF INTEGRATION OVER Z  — top-level worker functions
# ─────────────────────────────────────────────────────────────────────────────
def _init_intZ(shared_arrays, sp_list, noe, counter, total):
    global _arrZ, _species_listZ, _NOEqnZ, _counterZ, _totalZ
    _arrZ          = shared_arrays
    _species_listZ = sp_list
    _NOEqnZ        = noe
    _counterZ      = counter
    _totalZ        = total


def _worker_intZ(args):
    """Compute the Beta-PDF integral over Z for one (i, j) point."""
    i, j, a1, b1 = args

    import numpy as np
    from scipy.stats import beta as beta_dist

    with _counterZ.get_lock():
        _counterZ.value += 1
        done = _counterZ.value
    pct    = done / _totalZ * 100
    filled = int(40 * done / _totalZ)
    bar    = '█' * filled + '░' * (40 - filled)
    end    = '\n' if done == _totalZ else '\r'
    print(f"  Integration over Z   {done:>5}/{_totalZ}  {pct:5.1f}%  [{bar}]", end=end, flush=True)

    Z01_loc = _arrZ['Z01']
    CDF     = beta_dist.cdf(Z01_loc, a1, b1)

    def _intZ(f_last, df):
        return CDF[-1] * f_last - np.trapezoid(CDF[:, None, None] * df, Z01_loc, axis=0)

    result = {
        'i': i, 'j': j,
        'W':        _intZ(_arrZ['W1'][-1],        _arrZ['dWdZ']),
        'dH':       _intZ(_arrZ['dH1'][-1],       _arrZ['ddHdZ']),
        'T':        _intZ(_arrZ['T1'][-1],        _arrZ['dTdZ']),
        'Cp':       _intZ(_arrZ['Cp1'][-1],       _arrZ['dCpdZ']),
        'D':        _intZ(_arrZ['D1'][-1],        _arrZ['dDdZ']),
        'k':        _intZ(_arrZ['k1'][-1],        _arrZ['dkdZ']),
        'rho':      _intZ(_arrZ['rho1'][-1],      _arrZ['drhodZ']),
        'HRR':      _intZ(_arrZ['HRR1'][-1],      _arrZ['dHRRdZ']),
        'GYcYc':    _intZ(_arrZ['GYcYc1'][-1],    _arrZ['dGYcYcdZ']),
        'GYcZ':     _intZ(_arrZ['GYcZ1'][-1],     _arrZ['dGYcZdZ']),
        'GZYc':     _intZ(_arrZ['GZYc1'][-1],     _arrZ['dGZYcdZ']),
        'GZZ':      _intZ(_arrZ['GZZ1'][-1],      _arrZ['dGZZdZ']),
        'sourceYc': CDF[-1] * 0 - np.trapezoid(CDF[:, None, None] * _arrZ['dsourceYcdZ'],  Z01_loc, axis=0),
        'sourceYcV':CDF[-1] * 0 - np.trapezoid(CDF[:, None, None] * _arrZ['dsourceYcVdZ'], Z01_loc, axis=0),
        'sourceH2': CDF[-1] * 0 - np.trapezoid(CDF[:, None, None] * _arrZ['dsourceH2dZ'],  Z01_loc, axis=0),
    }
    if _NOEqnZ:
        result['sourceNO'] = _intZ(_arrZ['sourceNO1'][-1], _arrZ['dsourceNOdZ'])
    result['Yk'] = {sp: _intZ(_arrZ['Yk1'][sp][-1], _arrZ['dYkdZ'][sp]) for sp in _species_listZ}
    result['Dk'] = {sp: _intZ(_arrZ['Dk1'][sp][-1], _arrZ['dDkdZ'][sp]) for sp in _species_listZ}
    return result


if __name__ == "__main__":

    # Pack all arrays needed by C-integration workers into one dict
    shared_arrC = {
        'C': C,
        'W_fine': W_fine, 'dWdC': dWdC,
        'dH_fine': dH_fine, 'ddHdC': ddHdC,
        'D_fine': D_fine, 'dDdC': dDdC,
        'rho_fine': rho_fine, 'drhodC': drhodC,
        'k_fine': k_fine, 'dkdC': dkdC,
        'Cp_fine': Cp_fine, 'dCpdC': dCpdC,
        'HRR_fine': HRR_fine, 'dHRRdC': dHRRdC,
        'T_fine': T_fine, 'dTdC': dTdC,
        'GYcYc_fine': GYcYc_fine, 'dGYcYcdC': dGYcYcdC,
        'GYcZ_fine': GYcZ_fine, 'dGYcZdC': dGYcZdC,
        'GZYc_fine': GZYc_fine, 'dGZYcdC': dGZYcdC,
        'GZZ_fine': GZZ_fine, 'dGZZdC': dGZZdC,
        'sourceYc_fine': sourceYc_fine, 'dsourceYcdC': dsourceYcdC,
        'sourceYcV_fine': sourceYcV_fine, 'dsourceYcVdC': dsourceYcVdC,
        'sourceH2_fine': sourceH2_fine, 'dsourceH2dC': dsourceH2dC,
        'Yk_fine': Yk_fine, 'dYkdC': dYkdC,
        'Dk_fine': Dk_fine, 'dDkdC': dDkdC,
    }
    if NOEqn:
        shared_arrC['sourceNO_fine']  = sourceNO_fine
        shared_arrC['dsourceNOdC']    = dsourceNOdC

    # Build task list: one task per (i, j) pair
    tasks_C = [
        (i, j, float(a_C[i, j]), float(b_C[i, j]))
        for i in range(nCTilde - 2)
        for j in range(nCVar - 2)
    ]

    print(f"\nBeta-PDF integration over C — {len(tasks_C)} tasks on {N_WORKERS} workers...")

    counter_C = Value(ctypes.c_int, 0)

    WIntOverC        = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    DIntOverC        = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    dHIntOverC       = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    kIntOverC        = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    CpIntOverC       = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    rhoIntOverC      = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    TIntOverC        = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    HRRIntOverC      = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    GYcYcIntOverC    = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    GYcZIntOverC     = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    GZYcIntOverC     = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    GZZIntOverC      = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    sourceYcIntOverC = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    sourceYcVIntOverC= np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    sourceH2IntOverC = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    if NOEqn:
        sourceNOIntOverC = np.zeros((nZTilde, nCTilde - 2, nCVar - 2))
    YkIntOverC = {sp: np.zeros((nZTilde, nCTilde - 2, nCVar - 2)) for sp in species}
    DkIntOverC = {sp: np.zeros((nZTilde, nCTilde - 2, nCVar - 2)) for sp in species}

    resC_list = []
    with Pool(processes=N_WORKERS,
              initializer=_init_intC,
              initargs=(shared_arrC, species, NOEqn, counter_C, len(tasks_C))) as pool:
        for res in pool.imap_unordered(_worker_intC, tasks_C):
            resC_list.append(res)

    # Unpack results into pre-allocated arrays
    for res in resC_list:
        i, j = res['i'], res['j']
        WIntOverC[:, i, j]         = res['W']
        dHIntOverC[:, i, j]        = res['dH']
        DIntOverC[:, i, j]         = res['D']
        rhoIntOverC[:, i, j]       = res['rho']
        kIntOverC[:, i, j]    = res['k']
        CpIntOverC[:, i, j]        = res['Cp']
        HRRIntOverC[:, i, j]       = res['HRR']
        TIntOverC[:, i, j]         = res['T']
        GYcYcIntOverC[:, i, j]     = res['GYcYc']
        GYcZIntOverC[:, i, j]      = res['GYcZ']
        GZYcIntOverC[:, i, j]      = res['GZYc']
        GZZIntOverC[:, i, j]       = res['GZZ']
        sourceYcIntOverC[:, i, j]  = res['sourceYc']
        sourceYcVIntOverC[:, i, j] = res['sourceYcV']
        sourceH2IntOverC[:, i, j]  = res['sourceH2']
        if NOEqn:
            sourceNOIntOverC[:, i, j] = res['sourceNO']
        for sp in species:
            YkIntOverC[sp][:, i, j] = res['Yk'][sp]
            DkIntOverC[sp][:, i, j] = res['Dk'][sp]

    print("Integration over C completed!")

    # ── Add C-boundaries ──────────────────────────────────────────────────────
    WIntOverC01        = m.addBoundaries3D(WIntOverC,        WZTilde,        CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    dHIntOverC01       = m.addBoundaries3D(dHIntOverC,       dHZTilde,       CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    TIntOverC01        = m.addBoundaries3D(TIntOverC,        TZTilde,        CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    DIntOverC01        = m.addBoundaries3D(DIntOverC,        DZTilde,        CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    rhoIntOverC01      = m.addBoundaries3D(rhoIntOverC,      rhoZTilde,      CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    CpIntOverC01       = m.addBoundaries3D(CpIntOverC,       CpZTilde,       CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    kIntOverC01        = m.addBoundaries3D(kIntOverC,        kZTilde,        CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    HRRIntOverC01      = m.addBoundaries3D(HRRIntOverC,      HRRZTilde,      CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    GYcYcIntOverC01    = m.addBoundaries3D(GYcYcIntOverC,    GYcYcZTilde,    CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    GYcZIntOverC01     = m.addBoundaries3D(GYcZIntOverC,     GYcZZTilde,     CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    GZYcIntOverC01     = m.addBoundaries3D(GZYcIntOverC,     GZYcZTilde,     CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    GZZIntOverC01      = m.addBoundaries3D(GZZIntOverC,      GZZZTilde,      CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    sourceYcIntOverC01 = m.addBoundaries3D(sourceYcIntOverC, sourceYcZTilde, CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    sourceYcVIntOverC01= m.addBoundaries3D(sourceYcVIntOverC,sourceYcVZTilde,CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    sourceH2IntOverC01 = m.addBoundaries3D(sourceH2IntOverC, sourceH2ZTilde, CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    if NOEqn:
        sourceNOIntOverC01 = m.addBoundaries3D(sourceNOIntOverC, sourceNOZTilde, CZTilde, CTilde01, nZTilde, nCTilde, nCVar)
    YkIntOverC01 = {sp: m.addBoundaries3D(YkIntOverC[sp], YkZTilde[sp], CZTilde, CTilde01, nZTilde, nCTilde, nCVar) for sp in species}
    DkIntOverC01 = {sp: m.addBoundaries3D(DkIntOverC[sp], DkZTilde[sp], CZTilde, CTilde01, nZTilde, nCTilde, nCVar) for sp in species}

    print("Added boundaries for C = 0, C = 1, sigma_c = 0, sigma_c = sigma_c_MAX!")

    if plotTables:
        m.plot3DTable(CTILDE01, CVAR01, WIntOverC01[1, :, :],         " ", "W",                              r"\widetilde{C}", r"\widetilde{\sigma}_c")
        m.plot3DTable(CTILDE01, CVAR01, TIntOverC01[1, :, :],         " ", "T",                              r"\widetilde{C}", r"\widetilde{\sigma}_c")
        m.plot3DTable(CTILDE01, CVAR01, dHIntOverC01[1, :, :],        " ", r"\Delta\widetilde{h}_f^0",        r"\widetilde{C}", r"\widetilde{\sigma}_c")
        m.plot3DTable(CTILDE01, CVAR01, sourceYcIntOverC01[10, :, :], " ", r"\dot{\omega}_{H_2O}",            r"\widetilde{C}", r"\widetilde{\sigma}_c")
        m.plot3DTable(CTILDE01, CVAR01, sourceH2IntOverC01[1, :, :],  " ", r"\dot{\omega}_{H_2}",             r"\widetilde{C}", r"\widetilde{\sigma}_c")

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3 — Prepare for Z integration
    # ─────────────────────────────────────────────────────────────────────────
    nPOINTSZ = 1000
    Z01      = np.linspace(0, 1, nPOINTSZ)
    Z        = Z01[1:-1]

    WIntOverC01Reshaped        = WIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    TIntOverC01Reshaped        = TIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    DIntOverC01Reshaped        = DIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    dHIntOverC01Reshaped       = dHIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    kIntOverC01Reshaped        = kIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    CpIntOverC01Reshaped       = CpIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    rhoIntOverC01Reshaped      = rhoIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    HRRIntOverC01Reshaped      = HRRIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    GYcYcIntOverC01Reshaped    = GYcYcIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    GYcZIntOverC01Reshaped     = GYcZIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    GZYcIntOverC01Reshaped     = GZYcIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    GZZIntOverC01Reshaped      = GZZIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    sourceYcIntOverC01Reshaped = sourceYcIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    sourceYcVIntOverC01Reshaped= sourceYcVIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    sourceH2IntOverC01Reshaped = sourceH2IntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    if NOEqn:
        sourceNOIntOverC01Reshaped = sourceNOIntOverC01.reshape(nZTilde, nCTilde * nCVar).T
    YkIntOverC01Reshaped = {sp: YkIntOverC01[sp].reshape(nZTilde, nCTilde * nCVar).T for sp in species}
    DkIntOverC01Reshaped = {sp: DkIntOverC01[sp].reshape(nZTilde, nCTilde * nCVar).T for sp in species}

    # -- interpolate onto fine Z grid -----------------------------------------
    W2   = np.zeros((nCTilde * nCVar, nPOINTSZ))
    T2   = np.zeros((nCTilde * nCVar, nPOINTSZ))
    dH2  = np.zeros((nCTilde * nCVar, nPOINTSZ))
    D2   = np.zeros((nCTilde * nCVar, nPOINTSZ))
    k2   = np.zeros((nCTilde * nCVar, nPOINTSZ))
    Cp2  = np.zeros((nCTilde * nCVar, nPOINTSZ))
    rho2 = np.zeros((nCTilde * nCVar, nPOINTSZ))
    HRR2 = np.zeros((nCTilde * nCVar, nPOINTSZ))
    GYcYc2 = np.zeros((nCTilde * nCVar, nPOINTSZ))
    GYcZ2  = np.zeros((nCTilde * nCVar, nPOINTSZ))
    GZYc2  = np.zeros((nCTilde * nCVar, nPOINTSZ))
    GZZ2   = np.zeros((nCTilde * nCVar, nPOINTSZ))
    sourceYc2  = np.zeros((nCTilde * nCVar, nPOINTSZ))
    sourceYcV2 = np.zeros((nCTilde * nCVar, nPOINTSZ))
    sourceH22  = np.zeros((nCTilde * nCVar, nPOINTSZ))
    if NOEqn:
        sourceNO2 = np.zeros((nCTilde * nCVar, nPOINTSZ))
    Yk2 = {sp: np.zeros((nCTilde * nCVar, nPOINTSZ)) for sp in species}
    Dk2 = {sp: np.zeros((nCTilde * nCVar, nPOINTSZ)) for sp in species}

    for i in range(nCTilde * nCVar):
        W2[i, :]        = m.interpExt(Z01, ZTilde01, WIntOverC01Reshaped[i, :])
        T2[i, :]        = m.interpExt(Z01, ZTilde01, TIntOverC01Reshaped[i, :])
        dH2[i, :]       = m.interpExt(Z01, ZTilde01, dHIntOverC01Reshaped[i, :])
        D2[i, :]        = m.interpExt(Z01, ZTilde01, DIntOverC01Reshaped[i, :])
        k2[i, :]        = m.interpExt(Z01, ZTilde01, kIntOverC01Reshaped[i, :])
        Cp2[i, :]       = m.interpExt(Z01, ZTilde01, CpIntOverC01Reshaped[i, :])
        HRR2[i, :]      = m.interpExt(Z01, ZTilde01, HRRIntOverC01Reshaped[i, :])
        rho2[i, :]      = m.interpExt(Z01, ZTilde01, rhoIntOverC01Reshaped[i, :])
        sourceYc2[i, :] = m.interpExt(Z01, ZTilde01, sourceYcIntOverC01Reshaped[i, :])
        sourceYcV2[i,:] = m.interpExt(Z01, ZTilde01, sourceYcVIntOverC01Reshaped[i, :])
        sourceH22[i, :] = m.interpExt(Z01, ZTilde01, sourceH2IntOverC01Reshaped[i, :])
        if NOEqn:
            sourceNO2[i, :] = m.interpExt(Z01, ZTilde01, sourceNOIntOverC01Reshaped[i, :])
        GYcYc2[i, :] = m.interpExt(Z01, ZTilde01, GYcYcIntOverC01Reshaped[i, :])
        GYcZ2[i, :]  = m.interpExt(Z01, ZTilde01, GYcZIntOverC01Reshaped[i, :])
        GZYc2[i, :]  = m.interpExt(Z01, ZTilde01, GZYcIntOverC01Reshaped[i, :])
        GZZ2[i, :]   = m.interpExt(Z01, ZTilde01, GZZIntOverC01Reshaped[i, :])
        for sp in species:
            Yk2[sp][i, :] = m.interpExt(Z01, ZTilde01, YkIntOverC01Reshaped[sp][i, :])
            Dk2[sp][i, :] = m.interpExt(Z01, ZTilde01, DkIntOverC01Reshaped[sp][i, :])

    W1   = W2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    T1   = T2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    dH1  = dH2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    k1   = k2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    Cp1  = Cp2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    rho1 = rho2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    D1   = D2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    HRR1 = HRR2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    GYcYc1 = GYcYc2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    GYcZ1  = GYcZ2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    GZYc1  = GZYc2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    GZZ1   = GZZ2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    sourceYc1  = sourceYc2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    sourceYcV1 = sourceYcV2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    sourceH21  = sourceH22.T.reshape(nPOINTSZ, nCTilde, nCVar)
    if NOEqn:
        sourceNO1 = sourceNO2.T.reshape(nPOINTSZ, nCTilde, nCVar)
    Yk1 = {sp: Yk2[sp].T.reshape(nPOINTSZ, nCTilde, nCVar) for sp in species}
    Dk1 = {sp: Dk2[sp].T.reshape(nPOINTSZ, nCTilde, nCVar) for sp in species}

    # Z variance grid
    iexp      = np.arange(1, nZVar + 1, 1)
    ZVarMax01 = 0.25 * ((iexp - 1) / (nZVar - 1)) ** 2.5
    [ZVAR01, ZTILDE01] = np.meshgrid(ZVarMax01, ZTilde01)
    ZTILDE = ZTILDE01[1:-1, 1:-1]
    ZVAR   = ZVAR01[1:-1, 1:-1]
    ZVAR   = np.minimum(ZTILDE * (1 - ZTILDE), ZVAR)
    idx    = ZVAR == ZTILDE * (1 - ZTILDE)
    ZVAR[idx] = ZVAR[idx] - 1e-6

    ZZ_3d = np.tile(Z01[:, None, None], (1, nCTilde, nCVar))
    dZ    = np.diff(Z)

    a_Z = ZTILDE * (ZTILDE * (1 - ZTILDE) / ZVAR - 1)
    b_Z = (1 - ZTILDE) * (ZTILDE * (1 - ZTILDE) / ZVAR - 1)

    # Savitzky-Golay derivatives over Z
    dWdZ        = m.SavitzkyGolay(W1,        ZZ_3d, 0, 2, 1)
    ddHdZ       = m.SavitzkyGolay(dH1,       ZZ_3d, 0, 2, 1)
    dCpdZ       = m.SavitzkyGolay(Cp1,       ZZ_3d, 0, 2, 1)
    dTdZ        = m.SavitzkyGolay(T1,        ZZ_3d, 0, 2, 1)
    dDdZ        = m.SavitzkyGolay(D1,        ZZ_3d, 0, 2, 1)
    dkdZ        = m.SavitzkyGolay(k1,        ZZ_3d, 0, 2, 1)
    drhodZ      = m.SavitzkyGolay(rho1,      ZZ_3d, 0, 2, 1)
    dHRRdZ      = m.SavitzkyGolay(HRR1,      ZZ_3d, 0, 2, 1)
    dGYcYcdZ    = m.SavitzkyGolay(GYcYc1,    ZZ_3d, 0, 2, 1)
    dGYcZdZ     = m.SavitzkyGolay(GYcZ1,     ZZ_3d, 0, 2, 1)
    dGZYcdZ     = m.SavitzkyGolay(GZYc1,     ZZ_3d, 0, 2, 1)
    dGZZdZ      = m.SavitzkyGolay(GZZ1,      ZZ_3d, 0, 2, 1)
    dsourceYcdZ = m.SavitzkyGolay(sourceYc1, ZZ_3d, 0, 2, 1)
    dsourceYcVdZ= m.SavitzkyGolay(sourceYcV1,ZZ_3d, 0, 2, 1)
    dsourceH2dZ = m.SavitzkyGolay(sourceH21, ZZ_3d, 0, 2, 1)
    if NOEqn:
        dsourceNOdZ = m.SavitzkyGolay(sourceNO1, ZZ_3d, 0, 2, 1)
    dYkdZ = {sp: m.SavitzkyGolay(Yk1[sp], ZZ_3d, 0, 2, 1) for sp in species}
    dDkdZ = {sp: m.SavitzkyGolay(Dk1[sp], ZZ_3d, 0, 2, 1) for sp in species}

    # ─────────────────────────────────────────────────────────────────────────
    # PARALLEL BETA-PDF INTEGRATION OVER Z
    # ─────────────────────────────────────────────────────────────────────────
    shared_arrZ = {
        'Z01': Z01,
        'W1': W1, 'dWdZ': dWdZ,
        'dH1': dH1, 'ddHdZ': ddHdZ,
        'T1': T1, 'dTdZ': dTdZ,
        'Cp1': Cp1, 'dCpdZ': dCpdZ,
        'D1': D1, 'dDdZ': dDdZ,
        'k1': k1, 'dkdZ': dkdZ,
        'rho1': rho1, 'drhodZ': drhodZ,
        'HRR1': HRR1, 'dHRRdZ': dHRRdZ,
        'GYcYc1': GYcYc1, 'dGYcYcdZ': dGYcYcdZ,
        'GYcZ1': GYcZ1, 'dGYcZdZ': dGYcZdZ,
        'GZYc1': GZYc1, 'dGZYcdZ': dGZYcdZ,
        'GZZ1': GZZ1, 'dGZZdZ': dGZZdZ,
        'sourceYc1': sourceYc1, 'dsourceYcdZ': dsourceYcdZ,
        'sourceYcV1': sourceYcV1, 'dsourceYcVdZ': dsourceYcVdZ,
        'sourceH21': sourceH21, 'dsourceH2dZ': dsourceH2dZ,
        'Yk1': Yk1, 'dYkdZ': dYkdZ,
        'Dk1': Dk1, 'dDkdZ': dDkdZ,
    }
    if NOEqn:
        shared_arrZ['sourceNO1']   = sourceNO1
        shared_arrZ['dsourceNOdZ'] = dsourceNOdZ

    tasks_Z = [
        (i, j, float(a_Z[i, j]), float(b_Z[i, j]))
        for i in range(nZTilde - 2)
        for j in range(nZVar - 2)
    ]

    print(f"\nBeta-PDF integration over Z — {len(tasks_Z)} tasks on {N_WORKERS} workers...")

    counter_Z = Value(ctypes.c_int, 0)

    WIntOverZ        = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    dHIntOverZ       = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    TIntOverZ        = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    CpIntOverZ       = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    rhoIntOverZ      = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    DIntOverZ        = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    kIntOverZ        = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    HRRIntOverZ      = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    GYcYcIntOverZ    = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    GYcZIntOverZ     = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    GZYcIntOverZ     = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    GZZIntOverZ      = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    sourceYcIntOverZ = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    sourceYcVIntOverZ= np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    sourceH2IntOverZ = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    if NOEqn:
        sourceNOIntOverZ = np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar))
    YkIntOverZ = {sp: np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar)) for sp in species}
    DkIntOverZ = {sp: np.zeros((nZTilde - 2, nZVar - 2, nCTilde, nCVar)) for sp in species}

    resZ_list = []
    with Pool(processes=N_WORKERS,
              initializer=_init_intZ,
              initargs=(shared_arrZ, species, NOEqn, counter_Z, len(tasks_Z))) as pool:
        for res in pool.imap_unordered(_worker_intZ, tasks_Z):
            resZ_list.append(res)

    for res in resZ_list:
        i, j = res['i'], res['j']
        WIntOverZ[i, j]         = res['W']
        dHIntOverZ[i, j]        = res['dH']
        TIntOverZ[i, j]         = res['T']
        CpIntOverZ[i, j]        = res['Cp']
        DIntOverZ[i, j]         = res['D']
        kIntOverZ[i, j]         = res['k']
        rhoIntOverZ[i, j]       = res['rho']
        HRRIntOverZ[i, j]       = res['HRR']
        GYcYcIntOverZ[i, j]     = res['GYcYc']
        GYcZIntOverZ[i, j]      = res['GYcZ']
        GZYcIntOverZ[i, j]      = res['GZYc']
        GZZIntOverZ[i, j]       = res['GZZ']
        sourceYcIntOverZ[i, j]  = res['sourceYc']
        sourceYcVIntOverZ[i, j] = res['sourceYcV']
        sourceH2IntOverZ[i, j]  = res['sourceH2']
        if NOEqn:
            sourceNOIntOverZ[i, j] = res['sourceNO']
        for sp in species:
            YkIntOverZ[sp][i, j] = res['Yk'][sp]
            DkIntOverZ[sp][i, j] = res['Dk'][sp]

    # ── Add Z-boundaries ──────────────────────────────────────────────────────
    TIntOverZ01        = m.addBoundaries4D(TIntOverZ,        TIntOverC01,        TTable,              CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    dHIntOverZ01       = m.addBoundaries4D(dHIntOverZ,       dHIntOverC01,       dHTable,             CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    WIntOverZ01        = m.addBoundaries4D(WIntOverZ,        WIntOverC01,        MMWTable,            CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    DIntOverZ01        = m.addBoundaries4D(DIntOverZ,        DIntOverC01,        DTable,              CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    kIntOverZ01        = m.addBoundaries4D(kIntOverZ,   kIntOverC01,   kTable,         CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    CpIntOverZ01       = m.addBoundaries4D(CpIntOverZ,       CpIntOverC01,       CpTable,             CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    rhoIntOverZ01      = m.addBoundaries4D(rhoIntOverZ,      rhoIntOverC01,      rhoTable,            CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    HRRIntOverZ01      = m.addBoundaries4D(HRRIntOverZ,      HRRIntOverC01,      HRRTable,            CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    sourceYcIntOverZ01 = m.addBoundaries4D(sourceYcIntOverZ, sourceYcIntOverC01, omegakTable['H2O'],  CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    sourceYcVIntOverZ01= m.addBoundaries4D(sourceYcVIntOverZ,sourceYcVIntOverC01,sourceYcVTable,      CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    sourceH2IntOverZ01 = m.addBoundaries4D(sourceH2IntOverZ, sourceH2IntOverC01, omegakTable['H2'],   CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    if NOEqn:
        sourceNOIntOverZ01 = m.addBoundaries4D(sourceNOIntOverZ, sourceNOIntOverC01, omegakTable['NO'], CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    GYcYcIntOverZ01    = m.addBoundaries4D(GYcYcIntOverZ,    GYcYcIntOverC01,    GYcYcTable,          CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    GYcZIntOverZ01     = m.addBoundaries4D(GYcZIntOverZ,     GYcZIntOverC01,     GYcZTable,           CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    GZYcIntOverZ01     = m.addBoundaries4D(GZYcIntOverZ,     GZYcIntOverC01,     GZYcTable,           CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    GZZIntOverZ01      = m.addBoundaries4D(GZZIntOverZ,      GZZIntOverC01,      GZZTable,            CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar)
    YkIntOverZ01 = {sp: m.addBoundaries4D(YkIntOverZ[sp], YkIntOverC01[sp], YkTable[sp], CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar) for sp in species}
    DkIntOverZ01 = {sp: m.addBoundaries4D(DkIntOverZ[sp], DkIntOverC01[sp], DkTable[sp], CTable, CTilde01, nZTilde, nCTilde, nCVar, nZVar) for sp in species}
    
    print("Added boundaries for Z = 0, Z = 1, sigma_Z = 0, sigma_Z = sigma_Z_MAX!")

    # ── Permute  (Z, C, Cv, Zv) ───────────────────────────────────────────────
    TIntOverZ01T        = TIntOverZ01.transpose(0, 2, 3, 1)
    dHIntOverZ01T       = dHIntOverZ01.transpose(0, 2, 3, 1)
    DIntOverZ01T        = DIntOverZ01.transpose(0, 2, 3, 1)
    WIntOverZ01T        = WIntOverZ01.transpose(0, 2, 3, 1)
    CpIntOverZ01T       = CpIntOverZ01.transpose(0, 2, 3, 1)
    kIntOverZ01T        = kIntOverZ01.transpose(0, 2, 3, 1)
    rhoIntOverZ01T      = rhoIntOverZ01.transpose(0, 2, 3, 1)
    HRRIntOverZ01T      = HRRIntOverZ01.transpose(0, 2, 3, 1)
    sourceYcIntOverZ01T = sourceYcIntOverZ01.transpose(0, 2, 3, 1)
    sourceYcVIntOverZ01T= sourceYcVIntOverZ01.transpose(0, 2, 3, 1)
    sourceH2IntOverZ01T = sourceH2IntOverZ01.transpose(0, 2, 3, 1)
    if NOEqn:
        sourceNOIntOverZ01T = sourceNOIntOverZ01.transpose(0, 2, 3, 1)
    GYcYcIntOverZ01T    = GYcYcIntOverZ01.transpose(0, 2, 3, 1)
    GYcZIntOverZ01T     = GYcZIntOverZ01.transpose(0, 2, 3, 1)
    GZYcIntOverZ01T     = GZYcIntOverZ01.transpose(0, 2, 3, 1)
    GZZIntOverZ01T      = GZZIntOverZ01.transpose(0, 2, 3, 1)
    YkIntOverZ01T = {sp: YkIntOverZ01[sp].transpose(0, 2, 3, 1) for sp in species}
    DkIntOverZ01T = {sp: DkIntOverZ01[sp].transpose(0, 2, 3, 1) for sp in species}

    # ── Reshape for OpenFOAM ──────────────────────────────────────────────────
    def _reshape_oF(arr):
        return arr.reshape(nZTilde, nCTilde, -1, order='F').reshape(-1, 1)

    ToF         = _reshape_oF(TIntOverZ01T)
    dHoF        = _reshape_oF(dHIntOverZ01T)
    DoF         = _reshape_oF(DIntOverZ01T)
    WoF         = _reshape_oF(WIntOverZ01T)
    CpoF        = _reshape_oF(CpIntOverZ01T)
    koF         = _reshape_oF(kIntOverZ01T)
    rhooF       = _reshape_oF(rhoIntOverZ01T)
    HRRoF       = _reshape_oF(HRRIntOverZ01T)
    sourceYcoF  = _reshape_oF(sourceYcIntOverZ01T)
    sourceYcVoF = _reshape_oF(sourceYcVIntOverZ01T)
    sourceH2oF  = _reshape_oF(sourceH2IntOverZ01T)
    if NOEqn:
        sourceNOoF = _reshape_oF(sourceNOIntOverZ01T)
    GYcYcoF = _reshape_oF(GYcYcIntOverZ01T)
    GYcZoF  = _reshape_oF(GYcZIntOverZ01T)
    GZYcoF  = _reshape_oF(GZYcIntOverZ01T)
    GZZoF   = _reshape_oF(GZZIntOverZ01T)
    YkoF = {sp: _reshape_oF(YkIntOverZ01T[sp]) for sp in species}
    DkoF = {sp: _reshape_oF(DkIntOverZ01T[sp]) for sp in species}

    # ── Clip source terms ─────────────────────────────────────────────────────
    sourceYcoF[sourceYcoF < 0]   = 0
    sourceYcVoF[sourceYcoF < 0]  = 0
    if NOEqn:
        sourceNOoF[sourceNOoF < 0] = 0
    for sp in species:
        YkoF[sp][YkoF[sp] < 0] = 0

    ZVec    = ZTilde01
    CVec    = CTilde01
    ZVarVec = ZVarMax01
    CVarVec = CVarMax01

    # ── Export ────────────────────────────────────────────────────────────────
    if exportTable:
        print("Writing Table to file...\n")
        m.writeFGMTable4D(
            zLean, zRich, tin, p,
            Air["rhoAir"], Air["kAir"], Air["CpAir"], Air["WAir"], Air["DAir"], Air["D_NOAir"],
            Fuel["rhoFuel"], Fuel["kFuel"], Fuel["CpFuel"], Fuel["WFuel"], Fuel["DFuel"], Fuel["D_NOFuel"],
            ZTilde01, CTilde01, ZVarVec, CVarVec,
            sLVec.flatten(), tauVec.flatten(), lFVec.flatten(),
            sourceH2oF.flatten(), sourceNOoF.flatten() if NOEqn else np.zeros_like(sourceH2oF.flatten()),
            sourceYcoF.flatten(), sourceYcVoF.flatten(),
            ToF.flatten(), rhooF.flatten(), DoF.flatten(),
            DkoF['NO'].flatten(), koF.flatten(), CpoF.flatten(), WoF.flatten(),
            HRRoF.flatten(), dHoF.flatten(),
            YkoF['H'].flatten(), YkoF['H2'].flatten(), YkoF['H2O'].flatten(),
            YkoF['H2O2'].flatten(), YkoF['HO2'].flatten(),
            YkoF['O'].flatten(), YkoF['O2'].flatten(), YkoF['OH'].flatten(),
            YkoF['N2'].flatten(), YkoF['NO'].flatten(),
            GYcYcoF.flatten(), GYcZoF.flatten(), GZYcoF.flatten(), GZZoF.flatten(),
            nZTilde, nCTilde, nZVar, nCVar,
            nameMech, filename="fgmProperties"
        )
        print("Tables successfully written in --> fgmProperties\n")

    print("\nDone.")
