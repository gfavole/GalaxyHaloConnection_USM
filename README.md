# Connecting galaxies to their dark matter halos

**Practical session of the mini-course on the galaxy–halo connection**
Universidad Técnica Federico Santa María (USM), November 2026

Ginevra Favole — Instituto de Astrofísica de Canarias (IAC) and Universidad de La Laguna (ULL)
<gifavole@ull.edu.es>

---

This repository contains simple, commented implementations of the main empirical
models of the galaxy–halo connection. Each model populates the same AbacusSummit
dark-matter halos at z = 0.8 with galaxies and compares the clustering of the
resulting mocks with the **DESI DR1** measurements of luminous red galaxies (LRG3)
and emission-line galaxies (ELG1) at 0.8 < z < 1.1, and their cross-correlation.

Models 01–04 share the same workflow (`common/pipeline.py`):

```
DESI DR1 data  →  halo catalog  →  mock (the model)  →  FCFC: ξ0, ξ2, wp  →  χ² and figures
```

so that **the only thing that changes from one folder to the next is the model**.
Model 05, which needs galaxy properties beyond the halo, is tested against the
hydrodynamical simulation IllustrisTNG300 instead of DESI.

## The models

| Folder | Model | Main idea | Free parameters |
|---|---|---|---|
| [`01_HOD_standard`](01_HOD_standard) | HOD, Zheng et al. 2007 | Mean number of centrals and satellites as a function of halo mass. Complete sample: every massive halo has a central. | log M_min, σ_logM, log M_0, log M_1, α |
| [`02_HOD_Rocher23`](02_HOD_Rocher23) | Modified HOD for ELGs, Rocher et al. 2023 | Gaussian centrals, symmetric (GHOD) and asymmetric (mHMQ) (incomplete sample: ELGs avoid massive, quenched halos), velocity bias, strict conformity, modified satellite profile. | A_c, log M_c, σ_M, γ, A_s, log M_0, α, f_σv, f_sat |
| [`03_SHAM_standard`](03_SHAM_standard) | Standard SHAM on V_peak | Rank (sub)halos by V_peak and match them to a stellar mass function. No free occupation: only the scatter. | σ_AM |
| [`04_SHAM_modified`](04_SHAM_modified) | Modified / extended SHAM, Favole et al. 2016 | Incomplete samples: fixed satellite fraction, V_peak threshold, or a Gaussian probability in V_peak. | σ_AM, f_sat, V_peak,min or (V_mean, σ_V) |
| [`05_SHAM_decorated`](05_SHAM_decorated) | Decorated SHAM, Favole et al. 2022 | SHAM gives M_*; a secondary subhalo property (starvation redshift z_starve) sets the sSFR and colour at fixed M_*: assembly bias. Run on TNG300-1-Dark and compared with the hydrodynamical TNG300-1. | σ, secondary property, ρ |

Each folder has its own `README.md` with the equations, the parameters, the scripts
and some exercises.

## Installation

### 1. Python

Python ≥ 3.10 and a few standard packages:

```bash
git clone https://github.com/gfavole/GalaxyHaloConnection_USM.git
cd GalaxyHaloConnection_USM
pip install -r requirements.txt
```

### 2. FCFC

The correlation functions are computed with **FCFC** (C. Zhao,
<https://github.com/cheng-zhao/FCFC>), a fast C code for pair counting. It must be
compiled with OpenMP and, to read the mocks in HDF5, with HDF5 support:

```bash
git clone https://github.com/cheng-zhao/FCFC.git
cd FCFC
# in options.mk set:  WITH_HDF5 := T   (and the paths to HDF5 if needed)
make
```

This creates `FCFC_2PT_BOX` in the FCFC folder. On macOS, Apple's `clang` has no
OpenMP: install `gcc`, `libomp` and `hdf5` with Homebrew and set `CC` in `options.mk`
to the Homebrew `gcc` (for example `gcc-14`). If you cannot compile with HDF5, set
`CATALOG_FORMAT = "ascii"` in `config_common.py`: the mocks are then written as text.

### 3. Halo catalogs

Download the two course catalogs from the course Drive (link given in class) and
keep them on a disk with enough space:

| File | Content | Size |
|---|---|---|
| `course_halocenbasecat_z08.h5` | AbacusSummit c000 host halos at z = 0.8: positions, velocities, log M, R_vir, V_peak, V_max, concentration, ID | 7.4 GB |
| `course_halosatbasecat_z08.h5` | dark-matter particles of the halos (satellite candidates), with the ID of their host and a sampled V_peak | 2.1 GB |

For `05_SHAM_decorated`, download also the two TNG300 catalogs (a few tens of MB):
`TNG300-1-Dark_snap055_subhalos.h5` and `TNG300-1_snap055_galaxies.h5`.

The DESI DR1 measurements are already in the repository (`data/desi_dr1/`).

### 4. Configuration

```bash
cp config_common_example.py config_common.py
```

and edit the paths at the top of `config_common.py`: `HALO_CAT`, `PARTICLE_CAT` and
`FCFC_DIR` (models 01-04), and `TNG_HYDRO_CAT`, `TNG_DARK_CAT` (model 05). Every model
reads them from there; `config_common.py` is not tracked by git.

### 5. Check that everything works

Without the big catalogs, you can test the installation on small synthetic catalogs
(their clustering is not realistic):

```bash
python tools/make_test_catalogs.py test_catalogs/
# in config_common.py: HALO_CAT = "test_catalogs/test_halocen.h5",
#                      PARTICLE_CAT = "test_catalogs/test_halosat.h5", BOXSIZE = 300.
python 01_HOD_standard/run_hod.py
python -m pytest -q tests              # quick tests of the models, no FCFC needed
```

Remember to set the real paths and `BOXSIZE = 2000.` again afterwards.

## Running a model

Each folder has a `config.py` with the parameters of the model and one or more scripts:

```bash
cd 01_HOD_standard
python run_hod.py
```

The outputs go to `outputs/` inside the folder:

* `<model>_<T>.h5` — the mock (positions in redshift space along z, velocities, and
  model-specific columns such as log M_*, V_peak, log sSFR)
* `<model>_<T>_clustering.txt` — ξ0, ξ2 and wp of the mock
* `figures/<model>_<T>.pdf` — model vs DESI: s ξ0, s ξ2, rp wp
* `<model>_summary.json` — number of galaxies, satellite fraction and χ²

with T = `LRG3`, `ELG1` and `ELGxLRG` (cross-correlation).

On the full AbacusSummit box a run takes a few minutes on a laptop (reading the
catalogs takes most of the time) and needs a few GB of memory.

## Details common to all the models

* **Number of galaxies**: N = n(z_eff) · L³ · f_host, with n(z_eff = 0.95) from the DESI
  DR1 n(z) tables (ELG1: 6.57×10⁻⁴, LRG3: 2.87×10⁻⁴ h³ Mpc⁻³). The host catalog keeps a
  random fraction f_host = 0.4 of the halos (to make the files smaller); the mocks are
  diluted by the same fraction, so the clustering does not change.
* **Redshift space**: z_s = z + v_z / (aH) along the z axis (plane-parallel).
* **FCFC**: natural estimator in a periodic box, 26 log bins in s and r_p in
  [0.17, 32] h⁻¹ Mpc, 120 μ bins; wp with π_max = 40 h⁻¹ Mpc. The model is
  interpolated at the separations of the data.
* **χ²**: ξ0 and ξ2 with the full DESI covariance; wp with diagonal errors only.
  The range of scales is set by `SMIN`, `SMAX` in `config_common.py`.
* **Model band**: square root of the diagonal of mock covariances from small
  AbacusSummit boxes (500 h⁻¹ Mpc), with the Hartlap factor.
* **Colours**: LRG red, ELG blue, cross black.

## Repository structure

```
config_common_example.py   paths, simulation, clustering (copy to config_common.py)
common/                    catalogs, DESI data and χ², workflow, figures
clustering/                interface to FCFC and its .conf templates
hod/                       HOD models (Zheng, Rocher+23), population, quick fits
sham/                      SHAM: abundance matching, selections, PDF, decoration
data/desi_dr1/             DESI DR1 measurements, covariances, n(z)
01_HOD_standard/ ... 05_SHAM_decorated/    one folder per model
tests/                     quick tests (pytest)
tools/                     synthetic test catalogs (tools/tng/: fake TNG300 inputs of model 05)
```

## References

* Zheng Z. et al. 2007, ApJ 667, 760 (arXiv:astro-ph/0703457) — standard HOD
* Rocher A. et al. 2023, JCAP 10, 016 (arXiv:2306.06319) — HOD for DESI ELGs
* Kravtsov A. et al. 2004, ApJ 609, 35; Conroy C. et al. 2006, ApJ 647, 201 — SHAM
* Favole G. et al. 2016 (arXiv:1507.04356) — modified SHAM for ELGs
* Hearin A. & Watson D. 2013 (arXiv:1304.5557) — age matching
* Favole G. et al. 2022 (arXiv:2101.10733) — decorated SHAM in IllustrisTNG
* Nelson D. et al. 2019, Comput. Astrophys. Cosmol. 6, 2 — IllustrisTNG public data release
* Weaver J. et al. 2023, A&A 677, A184 — stellar mass function used in the SHAM
* Zhao C. — FCFC, <https://github.com/cheng-zhao/FCFC>
* DESI Collaboration 2025, DESI DR1 (arXiv:2503.14745)

If you use this code, please cite the papers listed above for the models, data and
tools you use.

© 2026 Ginevra Favole (IAC, ULL).
