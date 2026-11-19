# DESI DR1 — ELG1 and LRG3 (0.8 < z < 1.1)

Measurements used in Favole et al. 2026, ApJ 1002, 90, from the public DESI DR1 LSS catalogs.

| File | Content |
|---|---|
| `{T}_xi0_xi2.txt` | s, ξ0, ξ2, σ(ξ0), σ(ξ2) — 26 log bins in [0.17, 32] h⁻¹Mpc |
| `{T}_wp.txt` | r_p, w_p, σ(w_p) — π_max = 40 h⁻¹Mpc |
| `{T}_cov_mono.txt`, `{T}_cov_quad.txt` | data covariances of ξ0 and ξ2 (26×26); σ = √diag |
| `mock_covariances/{T}_cov_{mono,quad,wp}.txt` | covariances from mocks in small AbacusSummit boxes (500 h⁻¹Mpc); used only for the band around the model |
| `nz/{T}_{NGC,SGC}_nz.txt` | DESI DR1 n(z): zmid zlow zhigh n(z) Nbin Vol_bin |

T = `LRG3`, `ELG1`, `ELGxLRG` (cross-correlation). There is no data covariance for w_p:
its χ² uses the diagonal errors only.
