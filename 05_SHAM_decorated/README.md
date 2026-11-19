# 05 — Decorated SHAM on IllustrisTNG (Favole et al. 2022)

HOD and standard SHAM assume that galaxies depend only on halo mass (or V_peak). But
at fixed mass, halos that formed earlier are more clustered (**halo assembly bias**),
and if galaxy properties such as colour or star formation depend on how the halo grew,
galaxies inherit this **galaxy assembly bias**.

A model with secondary properties has to be tested where we know the right answer:
in a **hydrodynamical simulation**. Here the model runs on the dark-matter-only
simulation **TNG300-1-Dark** and is compared with the true galaxies of its
hydrodynamical twin **TNG300-1**, both at z = 0.82 in a 205 h⁻¹ Mpc box. This is the
"multi-population SHAM" of Favole et al. 2022 (arXiv:2101.10733), an extension of the
age matching of Hearin & Watson 2013 (arXiv:1304.5557).

## The model

1. **Standard SHAM**: the dark subhalos, ranked by V_peak, get the stellar masses of the
   TNG300-1 galaxies (the stellar mass function is that of the hydro), with a Gaussian
   scatter of 0.125 dex in M_* at fixed V_peak.
2. **Galaxy properties at fixed M_***: in bins of M_* (8.75–9.3, 9.3–10, 10–12.5 in
   log M_*/h⁻¹M☉), the conditional distributions P(log sSFR | M_*) and
   P(g − i | M_*) are measured in TNG300-1.
3. **Secondary matching**: each mock galaxy gets a value drawn from those distributions;
   the draws are sorted and assigned in order of a **secondary subhalo property**:
   the subhalos that stopped growing earliest get the lowest sSFR and the reddest colours.

   | Secondary property | Definition (SubLink main progenitor branch) |
   |---|---|
   | **z_starve** (default) | starvation redshift, max(z_1/2, z_char, z_acc): z_1/2 when the subhalo first had half of its peak mass, z_char when it first reached 10¹² h⁻¹M☉, z_acc when it became a subhalo for good |
   | c_infall | concentration at infall, V_max / V_vir when the subhalo was last a central |

   A parameter ρ controls the strength of the correlation (Gaussian copula,
   X′ = ρ z(X) + √(1−ρ²) N(0,1)): ρ = 1 is the full rank order of Favole+22,
   ρ = 0 removes the secondary dependence (no assembly bias).
4. **Samples**: the same selection for the hydro galaxies (truth) and for the mock, at the
   DESI DR1 densities: **LRG-like** = quenched (log sSFR < −10) with the largest M_*;
   **ELG-like** = star-forming (log sSFR > −9.5) with the largest SFR.

Note: in Favole+22 the conditional distributions are fitted with composite Gaussians;
here the mock values are drawn directly from the TNG300-1 galaxies of each M_* bin.

| Parameter (`config.py`) | Meaning |
|---|---|
| `sigma_logMstar` | scatter in M_* at fixed V_peak (0.125 dex) |
| `secondary` | `z_starve`, `c_infall` or `None` |
| `rho` | strength of the galaxy–halo correlation at fixed M_* |
| `mstar_bins` | M_* bins where the hydro distributions are measured |

## Input catalogs

Two small HDF5 files, given ready to use in the course Drive. They were built from the
public IllustrisTNG data (group catalogs and SubLink merger trees, Nelson et al. 2019):

| File | Content |
|---|---|
| `TNG300-1-Dark_snap055_subhalos.h5` | dark subhalos: positions, velocities, central/satellite, V_max, V_peak, z_1/2, z_char, z_acc, z_starve, c_infall |
| `TNG300-1_snap055_galaxies.h5` | hydro galaxies with log M_* > 8.75: positions, velocities, central/satellite, log M_*, SFR, log sSFR, g − i |

Set their paths (`TNG_DARK_CAT`, `TNG_HYDRO_CAT`) in `../config_common.py`. Box size, redshift and
Ω_m are read from the files. `tools/tng/make_fake_tng_inputs.py` writes fake files with
the same format to test the code.

## Scripts

| Script | What it does |
|---|---|
| `python run_decorated.py` | SHAM without secondary property, with z_starve and with c_infall, vs TNG300-1 |
| `python compare_rho.py` | z_starve with ρ = 0, 0.5, 1 vs TNG300-1 |

Outputs in `outputs/`: truth and model samples (`.h5`, with log M_*, log sSFR, g − i),
figures `figures/SHAM_decorated_vs_TNG_{LRG3,ELG1,ELGxLRG}.pdf` and a summary json with
the satellite fractions, median M_* and the large-scale ratio model / truth.

## Exercises

1. Run `run_decorated.py`. The SHAM without secondary property already has the right
   stellar masses and the right fraction of quenched galaxies at each M_*: why does it
   fail to reproduce the clustering of LRG-like and ELG-like galaxies?
2. Which secondary property, z_starve or c_infall, does better? For LRGs or for ELGs?
3. Run `compare_rho.py`: how much assembly bias is in TNG300?
4. Compare the satellite fractions of truth and model (summary json). Where do the
   differences come from? (Hint: subhalos lost in the dark-matter-only simulation, orphans.)
5. Try finer M_* bins (`mstar_bins`): does the model improve?
