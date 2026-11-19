# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Parameters of the standard HOD (Zheng et al. 2007, arXiv:astro-ph/0703457).
Paths, simulation and clustering settings come from ../config_common.py.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))           # the repository root: common/, hod/, clustering/
try:
    from config_common import *                # noqa: F401,F403  paths, box, cosmology, clustering
except ImportError:
    sys.exit("config_common.py not found: run  cp config_common_example.py config_common.py  "
             "in the repository root and edit the paths")

COLUMNS = HOD_COLUMNS                          # host halos: x, y, z, vx, vy, vz, logM, Rvir, host_id
PARTICLE_COLUMNS = HOD_PARTICLE_COLUMNS        # DM particles of each host
OUTDIR = str(HERE / "outputs")                 # mocks, measurements, figures

# ── Number density ──────────────────────────────────────────────────────────
# True: log M_min is adjusted at every run so that the mock has the DESI density.
TUNE_DENSITY = True

# ── Satellites ──────────────────────────────────────────────────────────────
# 'particles': on DM particles of their host (positions and velocities of the particles)
# 'nfw'      : on an analytic NFW profile with Gaussian velocities
# compare_placement.py runs both.
SAT_PLACEMENT = "particles"

# ── Zheng+07 parameters (masses in h^-1 Msun) ───────────────────────────────
#   log M_min   halo mass where <N_cen> = 1/2 (set by the density if TUNE_DENSITY)
#   sigma_logM  width of the central transition
#   log M_0     cut-off mass of the satellites
#   log M_1     mass with ~1 satellite
#   alpha       slope of the satellite power law
# Starting values to play with:
# ZHENG_LRG = dict(model="zheng", logMmin=12.9, sigma_logM=0.3, logM0=13.0, logM1=14.0, alpha=1.0)
# ZHENG_ELG = dict(model="zheng", logMmin=11.8, sigma_logM=0.5, logM0=11.5, logM1=13.0, alpha=0.8)
# Best fits to DESI DR1 with fit_hod.py (satellites on particles), shown in the lecture:
ZHENG_LRG = dict(model="zheng", logMmin=13.34, sigma_logM=0.7533, logM0=12.57, logM1=14.71, alpha=0.7559)
ZHENG_ELG = dict(model="zheng", logMmin=13.31, sigma_logM=0.9611, logM0=10.63, logM1=13.88, alpha=0.4969)

# Model used for each tracer
HOD = {"LRG3": ZHENG_LRG, "ELG1": ZHENG_ELG}
