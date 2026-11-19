# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
sham/decorated.py
=================
Decorated SHAM (multi-population SHAM) of Favole et al. 2022 (arXiv:2101.10733), an
extension of the age matching of Hearin & Watson 2013 (arXiv:1304.5557), applied to
IllustrisTNG: the model runs on the dark-matter-only simulation (TNG300-1-Dark) and is
compared with the galaxies of the hydrodynamical one (TNG300-1).

  1. Standard SHAM: the dark subhalos, ranked by V_peak, get the stellar masses of
     the hydrodynamical galaxies (same box: the stellar mass function is reproduced
     exactly), with a Gaussian scatter sigma_logMstar in M_* at fixed V_peak.
  2. In bins of M_*, the conditional distribution P(X_gal | M_*) of a galaxy property
     (log sSFR and g - i colour) is measured in the hydrodynamical simulation.
  3. Each mock galaxy gets a value drawn from P(X_gal | M_*); the draws are rank-ordered
     and assigned according to a SECONDARY SUBHALO PROPERTY X_halo (default z_starve):
     the subhalos that stopped growing earliest (highest z_starve) get the lowest sSFR /
     reddest colours. This keeps P(X_gal | M_*) of the hydro and adds the correlation
     between galaxy and halo at fixed M_*: galaxy assembly bias.
     The correlation strength is rho (Gaussian copula):
         X' = rho z(X_halo) + sqrt(1 - rho^2) N(0, 1)
     rho = 1: full rank order (Favole+22); rho = 0: properties independent of the halo
     (a colour-blind SHAM, no assembly bias).
  4. The same selection is applied to the hydro galaxies (truth) and to the mock:
     LRG-like = quenched (log sSFR < LOG_SSFR_Q) with the largest M_*;
     ELG-like = star-forming (log sSFR > LOG_SSFR_SF) with the largest SFR.
"""

import numpy as np
import pandas as pd
from scipy.special import ndtri

LOG_SSFR_Q, LOG_SSFR_SF = -10.0, -9.5       # yr^-1: quenched / star-forming thresholds


# ─────────────────────────────────────────────────────────────────────────────
# 1. SHAM with the stellar masses of the hydrodynamical simulation
# ─────────────────────────────────────────────────────────────────────────────
def sham_stellar_mass(vpeak, logmstar_hydro, sigma_logMstar, rng):
    """Assign the hydro stellar masses to subhalos ranked by V_peak.
    Without scatter: the k-th largest V_peak gets the k-th largest M_*.
    With scatter: log M_* + N(0, sigma) is re-ranked, so that the stellar mass function
    is still exactly that of the hydro, with ~sigma dex of scatter at fixed V_peak.
    Subhalos beyond the number of galaxies get NaN."""
    ms = np.sort(logmstar_hydro)[::-1]
    n = min(len(vpeak), len(ms))
    out = np.full(len(vpeak), np.nan)
    order = np.argsort(-vpeak)[:n]
    out[order] = ms[:n]                                     # no-scatter assignment
    if sigma_logMstar > 0:
        noisy = out[order] + rng.normal(0., sigma_logMstar, n)
        out[order[np.argsort(-noisy)]] = ms[:n]              # re-rank: same SMF, with scatter
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 2. Secondary property with a tunable correlation
# ─────────────────────────────────────────────────────────────────────────────
def gaussian_rank(x, rng):
    """Map x to a standard normal variable keeping its rank order (ties broken at random)."""
    order = np.lexsort((rng.uniform(size=len(x)), x))
    u = np.empty(len(x))
    u[order] = (np.arange(len(x)) + 0.5) / len(x)
    return ndtri(u)


def noisy_property(x, rho, rng):
    """X' = rho z(X) + sqrt(1 - rho^2) noise. NaN in x -> pure noise."""
    out = rng.normal(size=len(x))
    ok = np.isfinite(x)
    if ok.any() and rho != 0:
        out[ok] = rho * gaussian_rank(x[ok], rng) + np.sqrt(1. - rho**2) * out[ok]
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 3. Decoration: galaxy properties drawn from the hydro, ordered by the halo property
# ─────────────────────────────────────────────────────────────────────────────
def decorate(mock, hydro, secondary, rho, mstar_bins, rng, properties=('log_sSFR', 'g_i')):
    """Add the galaxy properties to the mock (DataFrame with logMstar and `secondary`).

    In each bin of M_*: draw as many values as mock galaxies from the hydro galaxies of
    that bin (P(X_gal | M_*)), sort them and assign them in order of X'(secondary):
    the highest X' (oldest halos) gets the lowest sSFR and the highest g - i."""
    mock = mock.copy()
    x = mock[secondary].values.astype(float) if secondary else np.full(len(mock), np.nan)
    xp = noisy_property(x, rho, rng)
    lm_m, lm_h = mock['logMstar'].values, hydro['logMstar'].values
    out = {prop: np.full(len(mock), np.nan) for prop in properties}
    for lo, hi in zip(mstar_bins[:-1], mstar_bins[1:]):
        im = np.flatnonzero((lm_m >= lo) & (lm_m < hi))
        ih = np.flatnonzero((lm_h >= lo) & (lm_h < hi))
        if len(im) == 0 or len(ih) == 0:
            continue
        draw = rng.choice(ih, size=len(im), replace=True)            # draws from P(X_gal | M_*)
        order = im[np.argsort(-xp[im])]                               # oldest halos first
        for prop in properties:
            vals = np.sort(hydro[prop].values[draw])                  # ascending
            # oldest halos get the reddest colour (largest g - i) and the lowest sSFR
            out[prop][order] = vals[::-1] if prop == 'g_i' else vals
    for prop in properties:
        mock[prop] = out[prop]
    return mock


# ─────────────────────────────────────────────────────────────────────────────
# 4. DESI-like selections, identical for the hydro galaxies and the mock
# ─────────────────────────────────────────────────────────────────────────────
def select_tracer(cat, tracer, n_target, log_ssfr_q=LOG_SSFR_Q, log_ssfr_sf=LOG_SSFR_SF):
    """LRG-like: the n_target quenched galaxies with the largest M_*.
    ELG-like: the n_target star-forming galaxies with the largest SFR (log sSFR + log M_*)."""
    ssfr = cat['log_sSFR'].values
    if tracer.startswith('LRG'):
        idx = np.flatnonzero(ssfr < log_ssfr_q)
        key = cat['logMstar'].values[idx]
    else:
        idx = np.flatnonzero(ssfr > log_ssfr_sf)
        key = ssfr[idx] + cat['logMstar'].values[idx]
    if len(idx) < n_target:
        raise ValueError(f"{tracer}: only {len(idx):,} candidates, need {n_target:,}")
    return cat.iloc[idx[np.argsort(-key)[:n_target]]].reset_index(drop=True)
