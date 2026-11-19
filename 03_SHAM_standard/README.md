# 03 — Standard SHAM on V_peak

Subhalo Abundance Matching assumes a **monotonic relation between a (sub)halo property
and a galaxy property**: the most massive galaxies live in the (sub)halos with the
deepest potential wells. It has no occupation function: centrals and satellites come
directly from host halos and subhalos.

## The model

1. **Scatter**: V′_peak = V_peak [1 + N(0, σ_AM)], where V_peak is the peak of the
   maximum circular velocity over the history of the (sub)halo (less affected by
   tidal stripping than the present mass).
2. **Stellar mass function**: double Schechter of Weaver et al. 2023 at 0.8 < z < 1.1,
   cumulative n(>M_*) = ∫ Φ d log M_*.
3. **Matching**: the (sub)halo with the k-th largest V′_peak gets the M_* with
   n(>M_*) = k / L³.
4. **Sample**: the N = n_DESI L³ objects with the largest M_* (or all with
   M_* ≥ `logMstar_min`).

Host halos give **centrals**, satellite candidates give **satellites**; each galaxy has
the position and velocity of its (sub)halo.

| Parameter | Meaning |
|---|---|
| σ_AM | scatter between V_peak and M_*: the **only** knob. Larger σ_AM lets lower-mass halos into the sample and lowers the large-scale bias |
| `logMstar_min` | optional stellar-mass threshold instead of the DESI density |

**About the satellite candidates of the course catalogs.** They are dark-matter
particles of each halo with a V_peak *sampled* from the subhalo V_peak function, not
real resolved subhalos. With pure ranking very few of them enter the sample
(f_sat ≈ 1%): see `04_SHAM_modified` for how to fix the satellite fraction.

**Key assumption**: the sample is *complete* in stellar mass. Does it hold for
LRGs? And for ELGs?

## Scripts

| Script | What it does |
|---|---|
| `python run_sham.py` | standard SHAM for LRG3, ELG1 and the cross |

The mocks `outputs/SHAM_standard_<T>.h5` also contain log M_* and V_peak.

## Exercises

1. Run with σ_AM = 0, 0.2, 0.35, 0.5: how do the large-scale bias and the satellite
   fraction change?
2. Why does the standard SHAM reproduce the LRG clustering better than the ELG one?
3. Compare the satellite fraction with that of the HOD fits of `01_HOD_standard`.
