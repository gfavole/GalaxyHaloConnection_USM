# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Parameters of the standard SHAM on V_peak. Paths and clustering settings come
from ../config_common.py.
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

COLUMNS = SHAM_COLUMNS                         # (sub)halos with V_peak
OUTDIR = str(HERE / "outputs")

# ── Standard SHAM ───────────────────────────────────────────────────────────
# sigma_AM      scatter on V_peak: V'_peak = V_peak [1 + N(0, sigma_AM)]. It is the ONLY
#               knob of the standard SHAM: larger scatter -> lower-mass halos enter the
#               sample -> lower large-scale bias.
# logMstar_min  None = each sample has the DESI number density (largest M_*);
#               a number = all galaxies with log M_* >= logMstar_min instead.
SHAM = {
    "LRG3": dict(model="standard", sigma_AM=0.35, logMstar_min=None),
    "ELG1": dict(model="standard", sigma_AM=0.35, logMstar_min=None),
}
