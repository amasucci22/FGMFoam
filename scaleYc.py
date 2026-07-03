# -*- coding: utf-8 -*-
"""
CREATE Yc Scaling function - Yc,b = f(Z) to scale C = Yc/Yc,b(Z)
Parallel version — each flamelet (phi) is solved in a separate worker process.
"""

import pandas as pd
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import os
import matplotlib as mpl
import ctypes
from multiprocessing import Pool, Value
import modules.functions as m

mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern']

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
save  = True
Soret = False

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
a         = 3           # controls curvature
s         = np.linspace(0, 1, nPointsZ)
phiLean   = 0.28
phiRich   = 5.0
phiVector = phiLean + (phiRich - phiLean) * (np.exp(a * s) - 1) / (np.exp(a) - 1)


# ─────────────────────────────────────────────────────────────────────────────
# SHARED COUNTER  (tracks progress across worker processes)
# ─────────────────────────────────────────────────────────────────────────────
def _init_counter(counter, total):
    global _counter, _total
    _counter = counter
    _total   = total


# ─────────────────────────────────────────────────────────────────────────────
# WORKER FUNCTION  (one call per phi, runs in a separate process)
# ─────────────────────────────────────────────────────────────────────────────
def runScaleYcFlamelet(args):
    """
    Solve one premixed flamelet and return the mixture fraction Z and the
    maximum progress-variable value YcMax at that phi.

    Must be a top-level function so multiprocessing can pickle it.
    All imports are done inside the function to avoid pickling issues with
    Cantera objects.
    """
    i, phi, p, tin, mech, Soret, cantera_dir = args

    import cantera as ct
    import numpy as np
    import modules.functions as m

    ct.add_directory(cantera_dir)

    # ── Solve flamelet ────────────────────────────────────────────────────────
    f, gas = m.premixedFlame(p, tin, phi, mech, Soret)

    with _counter.get_lock():
        _counter.value += 1
        done = _counter.value
    pct = done / _total * 100
    print(f"  Flamelets completed: {done}/{_total}  ({pct:.1f}%)", flush=True)

    # ── Fresh gas object at this phi ──────────────────────────────────────────
    gas2 = ct.Solution(mech)
    air  = "O2:0.21,N2:0.79"
    gas2.set_equivalence_ratio(phi=phi, fuel="H2:1", oxidizer=air)
    gas2.TPX = tin, p, gas2.X

    YcMax = float(np.max(f.Y[gas2.species_index('H2O')]))
    ZPhi = float(gas2.mixture_fraction(fuel="H2:1", oxidizer=air, element='Bilger'))

    return i, ZPhi, YcMax


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    ct.add_directory(CANTERA_DATA_DIR)

    p   = ct.one_atm
    tin = 298.15

    # ── Build argument list for workers ──────────────────────────────────────
    args_list = [
        (i, phi, p, tin, mech, Soret, CANTERA_DATA_DIR)
        for i, phi in enumerate(phiVector)
    ]

    # ── Run flamelets in parallel ─────────────────────────────────────────────
    print(f"\nRunning {nPointsZ} flamelets on {N_WORKERS} workers...\n")

    counter = Value(ctypes.c_int, 0)
    with Pool(processes=N_WORKERS,
              initializer=_init_counter,
              initargs=(counter, nPointsZ)) as pool:
        rawResults = pool.map(runScaleYcFlamelet, args_list)

    # Sort by original phi index to preserve order
    rawResults.sort(key=lambda x: x[0])

    # ── Unpack results ────────────────────────────────────────────────────────
    Z_values = np.array([r[1] for r in rawResults])
    YcMax    = np.array([r[2] for r in rawResults])
    YcMin    = np.zeros(nPointsZ)

    print("\nAll flamelets done. Building scaling table...")

    # ── Interpolate onto uniform Z grid ──────────────────────────────────────
    ZVec = np.linspace(0, 1, 1000)

    ZExt      = np.concatenate(([0], Z_values, [1]))
    YcMaxExt  = np.concatenate(([0], YcMax,    [0]))

    f          = interp1d(ZExt, YcMaxExt, fill_value='extrapolate')
    YcMaxInt   = f(ZVec)
    YcMinInt   = np.zeros_like(YcMaxInt)

    plt.plot(ZVec, YcMaxInt)
    plt.xlabel("Z")
    plt.ylabel("$Y_{c,b}$")
    plt.title("Scaling function $Y_{c,b}(Z)$")

    # ── Save ──────────────────────────────────────────────────────────────────
    df = pd.DataFrame({
        'Z':    ZVec,
        'YcMax': YcMaxInt,
        'YcMin': YcMinInt,
    })

    if save:
        os.makedirs(nameMech, exist_ok=True)
        df.to_csv(path + 'scaleYc' + nameMech + '.csv', index=False)
        m.scaleYc(ZVec, YcMinInt, YcMaxInt, filename="scalingYcTable")
        print("Scale Yc Function successfully saved in folder ---> " + nameMech + "!")