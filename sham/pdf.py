# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
sham/pdf.py
===========
Extended SHAM with a probability distribution in V_peak (Favole et al. 2016,
arXiv:1507.04356).

Idea: star-forming galaxies (ELGs) are an INCOMPLETE sample. They do not live in
the most massive (sub)halos, which host quenched galaxies, nor in the smallest
ones, which are too faint. Instead of keeping the N largest V_peak (standard SHAM),
each (sub)halo is selected with a probability that is a Gaussian in V'_peak:

    P(V'_peak) ∝ exp[ -(V'_peak - V_mean)^2 / (2 sigma_V^2) ]

and the satellite fraction f_sat is fixed: N_sat = f_sat N and N_cen = (1 - f_sat) N
are drawn separately, from the hosts and from the satellite candidates.

  V_mean   where the selection peaks: sets the typical halo mass, i.e. the
           LARGE-SCALE BIAS
  sigma_V  how broad the selection is (sigma_V -> infinity: random sample;
           sigma_V -> 0: a single V_peak)
  f_sat    sets the ONE-HALO TERM (small scales)

Optionally the same PDF can have different parameters for centrals and satellites
(V_mean_sat, sigma_V_sat); by default they share V_mean, sigma_V.
"""

import numpy as np
import pandas as pd

from sham.model import weighted_draw


def gaussian_pdf(v, v_mean, sigma_v):
    """Selection probability (unnormalised) of an object with V'_peak = v."""
    return np.exp(-0.5 * ((v - v_mean) / sigma_v)**2)


def select_pdf(halos, n_target, fsat, v_mean, sigma_v, rng, v_mean_sat=None, sigma_v_sat=None):
    """Draw n_target galaxies: (1 - fsat) n_target centrals among the hosts and
    fsat n_target satellites among the candidates, with probability P(V'_peak).
    halos must have Vpeakscatt (from abundance_matching) and is_sat."""
    if 'is_sat' not in halos:
        raise ValueError("the PDF selection needs an is_sat column")
    n_sat = int(round(fsat * n_target))
    pars = {0: (v_mean, sigma_v),
            1: (v_mean if v_mean_sat is None else v_mean_sat, sigma_v if sigma_v_sat is None else sigma_v_sat)}
    out = []
    for flag, n in ((0, int(n_target) - n_sat), (1, n_sat)):
        sub = halos[halos['is_sat'].values == flag]
        w = gaussian_pdf(sub['Vpeakscatt'].values, *pars[flag])
        out.append(sub.iloc[weighted_draw(w, n, rng)])
    return pd.concat(out, ignore_index=True)
