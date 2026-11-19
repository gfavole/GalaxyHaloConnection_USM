# 02 — Modified HOD for ELGs (Rocher et al. 2023)

ELGs are star-forming galaxies: they are an **incomplete** sample. Massive halos host
quenched centrals, so the probability of finding an ELG central does not saturate at
1, as in the standard HOD, but peaks and then drops. Rocher et al. 2023
(JCAP 10, 016; arXiv:2306.06319) fitted several such HODs to the DESI One-Percent
ELGs. Here the LRGs keep the standard Zheng HOD of `01_HOD_standard`.

## The model

**Centrals**: an asymmetric Gaussian in log M, the *modified High Mass Quenched*
model (mHMQ, their Eq. 3.4). The occupation peaks at M_c and has a tail towards high
masses controlled by the asymmetry γ:

$$\langle N_{\rm cen}\rangle = \frac{A_c}{\sqrt{2\pi}\,\sigma_M}\exp\left[-\frac{(\log M - \log M_c)^2}{2\sigma_M^2}\right]\left[1 + {\rm erf}\left(\frac{\gamma(\log M - \log M_c)}{\sqrt2\,\sigma_M}\right)\right]$$

With γ = 0 this is the **base model**, a symmetric Gaussian (GHOD, their Eq. 3.1).

**Satellites**, a power law independent of the centrals (log M_1 = 13 fixed):

$$\langle N_{\rm sat}\rangle = A_s\left(\frac{M - M_0}{M_1}\right)^{\alpha}$$

A_c and A_s are rescaled by a common factor to match the DESI DR1 density
(their mocks have n = 10⁻³ h³Mpc⁻³).

**Satellite physics** (always on the analytic NFW profile):

| Ingredient | Parameter | What it does |
|---|---|---|
| velocity bias | f_σv | satellite velocity dispersion = f_σv × that of the halo: changes ξ2 at small scales (Fingers of God) |
| strict conformity | `conformity=True` | a halo has ELG satellites only if it has an ELG central |
| modified profile | f_exp, τ, λ_NFW | a fraction f_exp follows dN/dr ∝ exp(−r/(τ r_s)) and can sit outside R_vir; the rest an NFW with r_s → r_s/λ_NFW |
| satellite fraction | f_sat | with conformity, A_c and A_s are set by the density **and** f_sat |

## Models in `config.py`

| Name | Centrals | Satellites | Values |
|---|---|---|---|
| `GHOD_VB` | symmetric Gaussian (base model) | velocity bias | Rocher+23 values, to be fitted to DESI DR1 |
| `MHMQ_VB` | asymmetric Gaussian (mHMQ) | velocity bias | best fit to DESI DR1 (lecture) |
| `MHMQ_VB_CONF` | asymmetric Gaussian (mHMQ) | velocity bias + strict conformity + modified profile (their best model) | best fit to DESI DR1 (lecture) |

`ELG_MODEL` chooses the model of `run_rocher.py`; `ELG_MODELS` the models compared by
`compare_elg_models.py` and `plot_hod_shapes.py`.

## Scripts

| Script | What it does |
|---|---|
| `python run_rocher.py` | LRG3 with Zheng, ELG1 with `ELG_MODEL`, and the cross |
| `python compare_elg_models.py` | ELG1 and ELG×LRG with Zheng, GHOD, mHMQ and mHMQ + conformity (same halos, same seed) |
| `python plot_hod_shapes.py` | ⟨N_cen⟩(M), ⟨N_sat⟩(M) of the models (no FCFC needed) |
| `python fit_hod.py ELG1 --model mhmq_vb` | quick fit (free: log M_c, σ_M, log M_0, α, f_σv; plus f_sat with conformity; γ fixed) |

Models of `fit_hod.py`: `ghod_vb`, `mhmq_vb`, `mhmq_vb_conf`, `zheng`.

## Exercises

1. Run `plot_hod_shapes.py`: how do the central occupations differ from Zheng?
   Where do the ELG satellites live in each model?
2. Run `compare_elg_models.py`: what improves most with respect to Zheng, the shape of
   the centrals or the satellite physics?
3. Set `f_sigma_v = 1` in `MHMQ_VB`: which statistic changes and at which scales?
4. Fit `ghod_vb` and compare it with `MHMQ_VB`: does the asymmetry γ matter for the clustering?
5. With strict conformity the same halo can host an LRG central and an ELG central
   (the two tracers are populated independently). Is this physical? What does it do
   to the ELG×LRG cross-correlation?
