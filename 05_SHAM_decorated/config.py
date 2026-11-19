# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Parameters of the decorated SHAM (Favole et al. 2022, arXiv:2101.10733) on IllustrisTNG:
the model runs on TNG300-1-Dark and is compared with the galaxies of TNG300-1 (hydro).
Paths (TNG catalogs, FCFC), DESI densities and clustering settings come from
../config_common.py.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
try:
    from config_common import *                # noqa: F401,F403
except ImportError:
    sys.exit("config_common.py not found: run  cp config_common_example.py config_common.py  "
             "in the repository root and edit the paths")
if 'TNG_HYDRO_CAT' not in globals():
    sys.exit("TNG_HYDRO_CAT / TNG_DARK_CAT missing in config_common.py: copy them from "
             "config_common_example.py and set the paths of the two TNG300 catalogs")

OUTDIR = str(HERE / "outputs")

# ── TNG input catalogs (course Drive) ──────────────────────────────────────
# Their paths, TNG_HYDRO_CAT and TNG_DARK_CAT, are set in ../config_common.py.
# The box size (205 h^-1 Mpc), redshift (0.82) and Omega_m are read from the files.
# Sample sizes: DESI DR1 densities n(z_eff) of config_common.py times the TNG volume.

# ── Decorated SHAM ──────────────────────────────────────────────────────────
# sigma_logMstar  scatter in log M_* at fixed V_peak (0.125 dex in Favole+22)
# secondary       subhalo property that orders the galaxy properties at fixed M_*:
#                 'z_starve' (starvation redshift), 'c_infall' (concentration at infall)
#                 or None (no secondary property)
# rho             strength of the correlation galaxy - halo property at fixed M_*
#                 (1 = full rank order as in Favole+22; 0 = no assembly bias)
# mstar_bins      bins of log M_* [h^-1 Msun] where P(sSFR | M_*) and P(g - i | M_*) are
#                 measured in the hydro (those of Favole+22; try finer bins)
DECORATED = dict(sigma_logMstar=0.125, secondary="z_starve", rho=1.0,
                 mstar_bins=(8.75, 9.3, 10.0, 12.5))

# run_decorated.py runs these models (same subhalos, same seed) and compares them with TNG300-1
MODELS = {
    "SHAM, no secondary property":  {**DECORATED, "secondary": None, "rho": 0.0},
    "decorated SHAM, z_starve":     DECORATED,
    "decorated SHAM, c_infall":     {**DECORATED, "secondary": "c_infall"},
}
