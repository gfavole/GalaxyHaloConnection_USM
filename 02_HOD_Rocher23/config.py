# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Parameters of the modified HOD for ELGs of Rocher et al. 2023 (JCAP 10, 016;
arXiv:2306.06319), fitted there to the DESI One-Percent ELGs. LRGs keep the
standard Zheng+07 HOD. Paths and clustering settings come from ../config_common.py.
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

COLUMNS = HOD_COLUMNS
PARTICLE_COLUMNS = HOD_PARTICLE_COLUMNS
OUTDIR = str(HERE / "outputs")
TUNE_DENSITY = True                # Zheng: log M_min; Rocher: amplitudes A_c, A_s
SAT_PLACEMENT = "particles"        # used by the Zheng LRGs; Rocher models force 'nfw'

# ── LRG: standard HOD (best fit to DESI DR1, as in 01_HOD_standard) ─────────
ZHENG_LRG = dict(model="zheng", logMmin=13.34, sigma_logM=0.7533, logM0=12.57, logM1=14.71, alpha=0.7559)
ZHENG_ELG = dict(model="zheng", logMmin=13.31, sigma_logM=0.9611, logM0=10.63, logM1=13.88, alpha=0.4969)

# ── ELG: Rocher+23 (masses in h^-1 Msun) ────────────────────────────────────
# Centrals: asymmetric Gaussian in log M, the modified High Mass Quenched model (mHMQ):
#   <N_cen> = A_c / (sqrt(2 pi) sigma_M) exp[-(log M - log M_c)^2 / (2 sigma_M^2)]
#             x [1 + erf(gamma (log M - log M_c) / (sqrt(2) sigma_M))]
# Satellites: A_s ((M - M_0) / M_1)^alpha, log M_1 = 13 fixed.
# A_c and A_s are rescaled to the DESI DR1 density (their mocks have n = 1e-3).
# Satellite physics (always on the analytic profile, placement='nfw'):
#   f_sigma_v               velocity bias
#   conformity              strict conformity: satellites only in halos with an ELG central
#   f_exp, tau, lambda_nfw  modified profile (a fraction f_exp on exp(-r / (tau r_s)))
#   fsat                    with conformity, A_c and A_s are set by the density AND f_sat
# Values: best fits to DESI DR1 with fit_hod.py, shown in the lecture
# (chi^2(xi0) ~ 128 and ~ 57, vs ~ 5e4 for Zheng).

# 1) mHMQ + velocity bias (their Table 3 model)
MHMQ_VB = dict(model="mhmq", Ac=0.03828, logMc=11.55, sigma_M=0.2737, gamma=7.06,
               As=0.003828, logM0=11.62, logM1=13.0, alpha=-0.1429,
               f_sigma_v=1.476, conformity=False, placement="nfw")

# 2) mHMQ + velocity bias + strict conformity + modified profile (their best model, Table 6)
MHMQ_VB_CONF = dict(model="mhmq", Ac=0.04193, logMc=11.39, sigma_M=0.6138, gamma=5.47,
                    As=0.2861, logM0=10.94, logM1=13.0, alpha=0.9758,
                    f_sigma_v=2.674, f_exp=0.58, tau=6.14, lambda_nfw=0.67,
                    conformity=True, fsat=0.01415, placement="nfw")

# 0) Base model: SYMMETRIC Gaussian (GHOD, their Eq. 3.1, i.e. gamma = 0) with the occupation
#    of their Table 2 and the velocity bias of their Table 3. Not fitted to DESI DR1:
#    fit it with  python fit_hod.py ELG1 --model ghod_vb
GHOD_VB = dict(model="ghod", Ac=1.0, logMc=11.89, sigma_M=0.11,
               As=0.08, logM0=11.72, logM1=13.0, alpha=-0.31,
               f_sigma_v=1.27, conformity=False, placement="nfw")

# ── What the scripts run ────────────────────────────────────────────────────
# run_rocher.py: LRG3 with Zheng, ELG1 with this model
ELG_MODEL = MHMQ_VB_CONF
HOD = {"LRG3": ZHENG_LRG, "ELG1": ELG_MODEL}
# compare_elg_models.py and plot_hod_shapes.py: these ELG1 models (same halos, same seed)
ELG_MODELS = {"Zheng+07": ZHENG_ELG,
              "Rocher+23 GHOD (symmetric) + v. bias": GHOD_VB,
              "Rocher+23 mHMQ + v. bias": MHMQ_VB,
              "Rocher+23 mHMQ + v. bias + conformity": MHMQ_VB_CONF}
