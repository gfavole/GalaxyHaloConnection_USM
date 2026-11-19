# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Parameters of the modified SHAM: fixed satellite fraction, V_peak threshold and the
Gaussian PDF in V_peak of Favole et al. 2016 (arXiv:1507.04356). Paths and
clustering settings come from ../config_common.py.
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

COLUMNS = SHAM_COLUMNS
OUTDIR = str(HERE / "outputs")

# All models first assign M_* by abundance matching on V'_peak = V_peak [1 + N(0, sigma_AM)];
# they differ in how the sample is selected.

# ── 1. Standard SHAM (reference) ────────────────────────────────────────────
SHAM_STANDARD = {
    "LRG3": dict(model="standard", sigma_AM=0.35),
    "ELG1": dict(model="standard", sigma_AM=0.35),
}

# ── 2. Modified SHAM: f_sat + V_peak threshold ──────────────────────────────
# fsat       fixed satellite fraction: centrals and satellite candidates are selected
#            separately (the pure ranking gives f_sat ~ 1% with the course catalogs)
# Vpeak_min  [km/s] random draw of the objects with V'_peak > Vpeak_min instead of the
#            largest V'_peak: an INCOMPLETE sample above a threshold. It sets the
#            large-scale bias.
# Vpeak_min_sat  [km/s] optional threshold for the satellites only. The satellite
#            candidates have much lower V_peak than the hosts (check_thresholds.py):
#            with one threshold for both, a high threshold leaves too few satellites
#            and a low one lowers the bias of the centrals.
# Values tuned by eye on the AbacusSummit course catalogs (Nov 2026): a single
# threshold of 350 km/s gave LRG3 a large-scale w_p ~30% below DESI, and f_sat = 0.10
# gave ELG1 a one-halo term ~3 times too high at r_p < 1 h^-1 Mpc.
SHAM_MODIFIED = {
    "LRG3": dict(model="modified", sigma_AM=0.35, fsat=0.04, Vpeak_min=450., Vpeak_min_sat=300.),
    "ELG1": dict(model="modified", sigma_AM=0.35, fsat=0.04, Vpeak_min=110.),
}

# ── 3. Extended SHAM with a Gaussian PDF in V_peak (Favole+16) ──────────────
# P(V'_peak) ~ exp[-(V'_peak - V_mean)^2 / (2 sigma_V^2)], separately for centrals and
# satellites with fixed f_sat. V_mean sets the large-scale bias, sigma_V the width of the
# selection, f_sat the one-halo term. STARTING VALUES: tune them against the data
# (exercise). Optional: V_mean_sat, sigma_V_sat for different satellite parameters.
SHAM_PDF = {
    # LRG3: with sigma_V = 80 the large-scale w_p was ~30% below DESI for V_mean = 550,
    # ~12% for 650 and ~7% for 720 (the many low-V_peak halos in the Gaussian tail
    # dominate the draw). A narrower sigma_V = 50 keeps fewer of them: ~6% below.
    "LRG3": dict(model="pdf", sigma_AM=0.35, fsat=0.04, V_mean=720., sigma_V=50.,
                 V_mean_sat=300., sigma_V_sat=80.),
    "ELG1": dict(model="pdf", sigma_AM=0.35, fsat=0.04, V_mean=150., sigma_V=40.),
}

# run_sham.py runs this one; compare_sham.py runs the three
SHAM = SHAM_PDF
