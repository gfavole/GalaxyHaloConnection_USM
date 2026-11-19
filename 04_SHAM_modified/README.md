# 04 — Modified and extended SHAM (Favole et al. 2016)

The standard SHAM keeps the N (sub)halos with the largest V_peak: a **complete** sample.
Star-forming galaxies are not: they avoid the most massive (sub)halos, where galaxies
are quenched. Favole et al. 2016 (arXiv:1507.04356) modified the SHAM for ELGs by
selecting (sub)halos with a probability in V_peak and fixing the satellite fraction.

All models first assign M_* by abundance matching on V′_peak = V_peak [1 + N(0, σ_AM)]
(as in `03_SHAM_standard`); they differ in **how the sample is selected**.

## The models in `config.py`

**1. `SHAM_STANDARD`** — the N largest V′_peak (reference).

**2. `SHAM_MODIFIED`** — fixed satellite fraction and V_peak threshold:

| Parameter | Meaning |
|---|---|
| `fsat` | N_sat = f_sat N satellites and (1 − f_sat) N centrals, selected separately: sets the **one-halo term** |
| `Vpeak_min` | random draw among the objects with V′_peak > V_peak,min instead of the largest V′_peak: an incomplete sample above a threshold that sets the **large-scale bias** (default 450 km/s for LRG, 110 km/s for ELG; check_thresholds.py lists the allowed values) |
| `Vpeak_min_sat` | optional threshold for the satellites only: the satellite candidates have much lower V_peak than the hosts, so one threshold for both leaves too few satellites or lowers the bias of the centrals (LRG3: 450 km/s for centrals, 300 km/s for satellites) |

**3. `SHAM_PDF`** — Gaussian probability distribution in V_peak (Favole+16):

$$P(V'_{\rm peak}) \propto \exp\left[-\frac{(V'_{\rm peak} - V_{\rm mean})^2}{2\sigma_V^2}\right]$$

applied separately to centrals and satellites with fixed f_sat.

| Parameter | Meaning |
|---|---|
| V_mean | where the selection peaks: typical halo mass → **large-scale bias** |
| σ_V | width of the selection (σ_V → ∞: random sample) |
| f_sat | **one-halo term** |
| `V_mean_sat`, `sigma_V_sat` | optional, different PDF for the satellites |

The values in `config.py` are **starting values**: tuning them against DESI is the
exercise.

## Scripts

| Script | What it does |
|---|---|
| `python check_thresholds.py` | how many centrals and satellite candidates have V'_peak above each threshold, and the largest f_sat each threshold allows (no FCFC) |
| `python run_sham.py` | runs `config.SHAM` (by default `SHAM_PDF`) for LRG3, ELG1 and the cross |
| `python compare_sham.py` | standard vs modified vs PDF, same (sub)halos and seed |

## Exercises

1. Run `compare_sham.py`. Which statistic and which scales does each modification change?
2. With `SHAM_MODIFIED`, change `fsat` keeping `Vpeak_min` fixed, then the opposite:
   check that f_sat moves the small scales and V_peak,min the large scales.
3. Tune V_mean, σ_V and f_sat of the ELGs to reproduce ELG1. What halo masses do the
   ELGs occupy in your best model? (The mock contains V_peak: compare with the HOD
   of `02_HOD_Rocher23`.)
4. Why is a Gaussian in V_peak similar to the Gaussian central HOD of Rocher+23?
