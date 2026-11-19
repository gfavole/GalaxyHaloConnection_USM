# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Before choosing f_sat and V_peak_min for the modified SHAM: how many centrals
(host halos) and satellite candidates have V'_peak above each threshold, after the
same scatter sigma_AM as the model. No FCFC needed.

The modified SHAM needs, for each tracer,
    (1 - f_sat) N_target  centrals   with V'_peak > V_peak_min
         f_sat  N_target  satellites with V'_peak > V_peak_min
so the largest f_sat allowed by a threshold is N_sat(> V_peak_min) / N_target.

    python check_thresholds.py
"""
import numpy as np

import config
from common.catalogs import read_halos
from common.pipeline import n_target_of
from sham.mock import objects_of

THRESHOLDS = (80, 100, 110, 130, 150, 200, 250, 300, 350, 370, 400, 450, 500)   # km/s

if __name__ == '__main__':
    halos = read_halos(config.HALO_CAT, config.COLUMNS, boxsize=config.BOXSIZE)
    objs = objects_of(halos, config)
    rng = np.random.default_rng(getattr(config, 'SEED', 42))
    for T, p in config.SHAM_MODIFIED.items():
        _, n_target = n_target_of(config, T)
        sigma = p.get('sigma_AM', 0.)
        vs = objs['Vpeak'].values * (1. + rng.normal(0., abs(sigma), len(objs)))   # V'_peak
        sat = objs['is_sat'].values == 1
        vc, vsat = np.sort(vs[~sat]), np.sort(vs[sat])
        print(f"\n{T}: N_target = {n_target:,}   (sigma_AM = {sigma};"
              f" now f_sat = {p.get('fsat')}, V_peak_min = {p.get('Vpeak_min')} km/s)")
        print(f"   {'V_min':>6} {'N_cen(>V)':>12} {'N_sat(>V)':>12} {'max f_sat':>10}  centrals enough?")
        for v in THRESHOLDS:
            nc = len(vc) - np.searchsorted(vc, v)
            ns = len(vsat) - np.searchsorted(vsat, v)
            print(f"   {v:>6} {nc:>12,} {ns:>12,} {min(ns / n_target, 1.):>10.3f}  "
                  f"{'yes' if nc >= n_target else 'no'}")
