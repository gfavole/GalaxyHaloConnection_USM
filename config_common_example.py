# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Settings shared by ALL the models of the course: paths, simulation, clustering.

    cp config_common_example.py config_common.py      # then edit the three paths below

Each model folder has its own config.py with the model parameters; it imports
everything defined here.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ── Paths: EDIT THESE ────────────────────────────────────────────────────────
# Course catalogs (AbacusSummit c000, z = 0.8), downloaded from the course Drive
HALO_CAT     = "/path/to/course_halocenbasecat_z08.h5"    # host halos (centrals)
PARTICLE_CAT = "/path/to/course_halosatbasecat_z08.h5"    # DM particles = satellite candidates
FCFC_DIR     = "/path/to/FCFC/"                           # folder containing FCFC_2PT_BOX
# IllustrisTNG300 catalogs at z = 0.82, for 05_SHAM_decorated only (course Drive)
TNG_HYDRO_CAT = "/path/to/TNG300-1_snap055_galaxies.h5"   # hydro galaxies: M_*, sSFR, g - i
TNG_DARK_CAT  = "/path/to/TNG300-1-Dark_snap055_subhalos.h5"  # dark subhalos: V_peak, z_starve, c_infall

DATA_DIR = str(ROOT / "data/desi_dr1")                    # DESI DR1 measurements (in the repo)

# ── Catalog columns: internal name -> dataset name in the HDF5 files ────────
# HOD: host halos with mass, radius and ID (the ID links the particles to their host)
HOD_COLUMNS = dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz",
                   logM="logMcen", Rvir="Rvir", host_id="IDcen")
HOD_PARTICLE_COLUMNS = dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz", host_id="IDcen")
RVIR_TO_MPC = 1.                     # Rvir is already in h^-1 Mpc in the course catalogs
# SHAM: (sub)halos with V_peak (the decorated SHAM reads its own TNG catalogs, see 05_SHAM_decorated)
SHAM_COLUMNS = dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz", Vpeak="Vpeak")
SHAM_PARTICLE_COLUMNS = dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz", Vpeak="Vpeak")

# ── Simulation and cosmology (AbacusSummit c000) ─────────────────────────────
BOXSIZE = 2000.                      # h^-1 Mpc
OM      = 0.31519
Z_EFF   = 0.95                       # redshift used for n(z_eff) of DESI and for the RSD
HOST_FRACTION = 0.4                  # random fraction of halos kept in the host catalog
                                     # (overridden by the 'host_fraction' attribute of the file);
                                     # mocks have N = n_DESI L^3 HOST_FRACTION galaxies

# ── Clustering ───────────────────────────────────────────────────────────────
TRACERS = ("LRG3", "ELG1")           # DESI DR1 samples at 0.8 < z < 1.1
CROSS   = True                       # also the ELG x LRG cross-correlation
CATALOG_FORMAT = "hdf5"              # "hdf5" (FCFC compiled with WITH_HDF5 := T) or "ascii"
FCFC_NTHREADS  = 4                   # OpenMP threads for FCFC
SMIN, SMAX = None, None              # range of s (r_p) used in the chi^2 [h^-1 Mpc]; None = all
SEED = 42                            # random seed: same seed, same mock
