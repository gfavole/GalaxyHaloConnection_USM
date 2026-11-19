# 01 — Standard HOD (Zheng et al. 2007)

The Halo Occupation Distribution describes the **mean number of galaxies per halo as
a function of halo mass only**, separately for centrals and satellites.

## The model

$$\langle N_{\rm cen}\rangle(M) = \frac12\left[1 + {\rm erf}\left(\frac{\log M - \log M_{\rm min}}{\sigma_{\log M}}\right)\right]$$

$$\langle N_{\rm sat}\rangle(M) = \langle N_{\rm cen}\rangle(M)\left(\frac{M - M_0}{M_1}\right)^{\alpha} \quad (M > M_0)$$

| Parameter | Meaning | Main effect on the clustering |
|---|---|---|
| log M_min | mass where half of the halos have a central | large-scale bias (here set by the DESI density) |
| σ_logM | width of the central transition | how many low-mass halos enter: bias |
| log M_0 | minimum mass of halos with satellites | small scales |
| log M_1 | mass of a halo with ~1 satellite | satellite fraction: one-halo term |
| α | slope of the satellite occupation | which halos host the satellites |

* **Centrals**: one Bernoulli draw per halo with probability ⟨N_cen⟩, at the halo
  centre with the halo velocity.
* **Satellites**: Poisson number with mean ⟨N_sat⟩. Two ways of placing them
  (`SAT_PLACEMENT` in `config.py`):
  * `particles` — on dark-matter particles of the same halo (linked by the host ID),
    with the particle positions and velocities;
  * `nfw` — on an NFW profile truncated at R_vir (concentration of Duffy et al. 2008)
    with the halo velocity plus an isotropic Gaussian dispersion σ₁D = V_vir/√2.
* **Density**: with `TUNE_DENSITY = True`, log M_min is adjusted so that the mock has
  exactly the DESI number density; the other parameters are set in `config.py`.

**Key assumption**: the sample is *complete*: ⟨N_cen⟩ → 1 in massive halos.
This works for LRGs (massive, red, mass-selected galaxies), but not for ELGs
(star-forming galaxies, absent from the most massive halos): see `02_HOD_Rocher23`.

## Scripts

| Script | What it does |
|---|---|
| `python run_hod.py` | runs the HOD in `config.py` for LRG3, ELG1 and the cross |
| `python compare_placement.py` | same HOD, satellites on particles vs NFW (same seed) |
| `python fit_hod.py LRG3` | quick Nelder-Mead fit to DESI (free: σ_logM, log M_0, log M_1, α) |

Options of `fit_hod.py`: `--smin 5` (only s ≥ 5 h⁻¹Mpc), `--stats xi0 wp`, `--maxeval 30`.
The best fit is printed in `config.py` format and saved in `outputs/fit/`.

`config.py` contains the best fits to DESI DR1 shown in the lecture
(satellites on particles) and, commented, simpler starting values.

## Exercises

1. Change one parameter at a time (σ_logM, log M_1, α) and look at ξ0, ξ2 and wp:
   which scales does each parameter move?
2. Run `compare_placement.py`. Why is ξ2 at small scales so different between
   particles and NFW? (Hint: velocities, Fingers of God.)
3. Fit LRG3 and ELG1 with `fit_hod.py`. Why does the fit work for LRGs and not for
   ELGs? Look at the best-fit σ_logM and log M_min of the ELGs.
4. Fit only large scales (`--smin 5`): does the ELG fit improve? What does it tell
   you about the one-halo term?
