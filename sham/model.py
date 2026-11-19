# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
sham/model.py
=============
Subhalo Abundance Matching (SHAM) on V_peak, the peak of the maximum circular
velocity over the whole history of each (sub)halo.

Standard SHAM (Kravtsov et al. 2004; Conroy et al. 2006; Trujillo-Gomez et al. 2011):

  1. Scatter:   V'_peak = V_peak [1 + N(0, sigma_AM)]
  2. Galaxies:  cumulative number density from a stellar mass function,
                n_gal(>M_*) = int_{M_*} Phi(M) dlog M
  3. Matching:  the (sub)halo with the k-th largest V'_peak gets the M_* with
                n_gal(>M_*) = k / V_box          (rank order, one galaxy per (sub)halo)
  4. Sample:    the N_target objects with the largest M_* (density of the data),
                or all the objects above a stellar-mass threshold.

Hosts give centrals and subhalos (here: satellite candidates) give satellites;
each galaxy keeps the position and velocity of its (sub)halo. The model assumes a
COMPLETE, mass-selected sample: the only free parameter is sigma_AM.

Modified SHAM (selection of the sample, see select()):
  fsat        the satellite fraction is fixed: centrals and satellites are chosen
              separately (Favole et al. 2016, arXiv:1507.04356)
  vpeak_min   INCOMPLETE sample above a threshold: a random draw of the objects
              with V'_peak > vpeak_min instead of the largest V'_peak
The Gaussian PDF in V_peak of Favole+16 is in sham/pdf.py.
"""

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid


# ─────────────────────────────────────────────────────────────────────────────
# 1. Stellar mass function
# ─────────────────────────────────────────────────────────────────────────────
def smf_weaver23(logM):
    """Double Schechter SMF of Weaver et al. 2023 (COSMOS2020, 0.8 < z < 1.1),
    Phi in Mpc^-3 dex^-1, M_* in Msun (h = 0.7 convention of the paper)."""
    phi1, phi2 = 0.70e-3, 0.77e-3
    logmstar, a1, a2 = 11.01, -1.35, -0.72
    x = 10**(logM - logmstar)
    return np.log(10) * np.exp(-x) * (phi1 * x**(1 + a1) + phi2 * x**(1 + a2))


def cumulative_smf(smf=smf_weaver23, logM_min=8., logM_max=13., h=0.7):
    """(log M grid, n(>log M) in h^3 Mpc^-3)."""
    lm = np.linspace(logM_min, logM_max, 2001)
    phi = smf(lm) / h**3                               # Mpc^-3 -> h^3 Mpc^-3
    ncum = cumulative_trapezoid(phi[::-1], lm[::-1], initial=0)[::-1] * -1
    return lm, ncum


# ─────────────────────────────────────────────────────────────────────────────
# 2. Abundance matching
# ─────────────────────────────────────────────────────────────────────────────
def abundance_matching(halos, sigma_AM, boxsize, rng, smf=smf_weaver23, h_smf=0.7):
    """Add the columns Vpeakscatt (V'_peak) and logMstar (NaN below the SMF range)."""
    halos = halos.copy()
    # step 1: scatter on V_peak
    vs = halos['Vpeak'].values * (1. + rng.normal(0., abs(sigma_AM), len(halos)))
    # step 2-3: rank in V'_peak and match to the cumulative SMF
    order = np.argsort(vs)[::-1]
    n_rank = np.arange(1, len(halos) + 1) / boxsize**3
    lm, ncum = cumulative_smf(smf, h=h_smf)
    ok = n_rank <= ncum[0]
    logMstar = np.full(len(halos), np.nan)
    # ncum decreases with log M: interpolate log M(n) on the reversed arrays
    logMstar[order[ok]] = np.interp(n_rank[ok], ncum[::-1], lm[::-1])
    halos['Vpeakscatt'] = vs
    halos['logMstar'] = logMstar
    return halos


# ─────────────────────────────────────────────────────────────────────────────
# 3. Selection of the galaxy sample
# ─────────────────────────────────────────────────────────────────────────────
def weighted_draw(weights, n, rng):
    """Indices of n objects drawn WITHOUT repetition with probability proportional to
    `weights` (Efraimidis & Spirakis 2006: the n smallest -log(u) / w). Fast for
    tens of millions of objects."""
    w = np.asarray(weights, dtype=float)
    good = np.flatnonzero(w > 0)
    if len(good) < n:
        raise ValueError(f"only {len(good):,} objects with non-zero probability, need {n:,}. "
                         "With a V_peak threshold: lower V_peak_min or f_sat "
                         "(04_SHAM_modified/check_thresholds.py lists the allowed values)")
    # log of the key -log(u) / w (in log space to avoid overflows for tiny weights)
    logkey = np.log(-np.log(rng.uniform(size=len(good)))) - np.log(w[good])
    return good[np.argpartition(logkey, n - 1)[:n]]


def _pick(sub, n, vpeak_min, rng):
    """n objects of sub: the largest V'_peak (standard), or a random draw among those
    with V'_peak > vpeak_min (incomplete sample above a threshold)."""
    if vpeak_min is None:
        idx = np.argsort(sub['Vpeakscatt'].values)[::-1][:n]
        return sub.iloc[idx]
    w = (sub['Vpeakscatt'].values > vpeak_min).astype(float)        # flat PDF above the threshold
    return sub.iloc[weighted_draw(w, n, rng)]


def select(halos, n_target=None, logMstar_min=None, fsat=None, vpeak_min=None, rng=None,
           vpeak_min_sat=None):
    """Galaxy sample from the abundance-matched objects.

    Standard SHAM : the n_target objects with the largest M_* (largest V'_peak),
                    or all objects with M_* >= logMstar_min.
    fsat          : fsat n_target satellites and (1 - fsat) n_target centrals,
                    chosen separately (needs the is_sat column).
    vpeak_min     : a random draw of the objects with V'_peak > vpeak_min instead
                    of the largest V'_peak (this sets the large-scale bias).
    vpeak_min_sat : threshold for the satellites only (default: vpeak_min). The
                    satellite candidates have lower V_peak than the hosts, so a
                    single threshold can leave too few satellites.
    """
    if (n_target is None) == (logMstar_min is None):
        raise ValueError("give exactly one of n_target, logMstar_min")
    if logMstar_min is not None:
        return halos[halos['logMstar'] >= logMstar_min].reset_index(drop=True)
    rng = rng if rng is not None else np.random.default_rng()
    if fsat is None:
        return _pick(halos, int(n_target), vpeak_min, rng).reset_index(drop=True)
    if 'is_sat' not in halos:
        raise ValueError("fsat needs an is_sat column")
    n_sat = int(round(fsat * n_target))
    v_sat = vpeak_min if vpeak_min_sat is None else vpeak_min_sat
    out = [_pick(halos[halos['is_sat'].values == flag], n, v, rng)
           for flag, n, v in ((0, int(n_target) - n_sat, vpeak_min), (1, n_sat, v_sat))]
    return pd.concat(out, ignore_index=True)
