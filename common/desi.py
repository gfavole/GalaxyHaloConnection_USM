# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
common/desi.py
==============
DESI DR1 measurements, number densities, model errors and chi^2.

Files in data/desi_dr1/ (see its README), with T = LRG3, ELG1 or ELGxLRG:
  {T}_xi0_xi2.txt     s  xi0  xi2  err(xi0)  err(xi2)
  {T}_wp.txt          rp wp  err(wp)
  {T}_cov_mono.txt, {T}_cov_quad.txt            data covariances of xi0 and xi2
  mock_covariances/{T}_cov_{mono,quad,wp}.txt   mock covariances (band around the model)
  nz/{T}_{NGC,SGC}_nz.txt                       DESI DR1 n(z)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate


# ─────────────────────────────────────────────────────────────────────────────
# 1. Clustering data
# ─────────────────────────────────────────────────────────────────────────────
def load_clustering(data_dir, tracer):
    """DESI DR1 xi0(s), xi2(s) and wp(rp) with their diagonal errors."""
    data_dir = Path(data_dir)
    x = np.loadtxt(data_dir / f"{tracer}_xi0_xi2.txt", comments="#")
    w = np.loadtxt(data_dir / f"{tracer}_wp.txt", comments="#")
    return dict(s=x[:, 0], xi0=x[:, 1], xi2=x[:, 2], e0=x[:, 3], e2=x[:, 4],
                rp=w[:, 0], wp=w[:, 1], ewp=w[:, 2])


# ─────────────────────────────────────────────────────────────────────────────
# 2. Number density: it sets how many galaxies each mock has
# ─────────────────────────────────────────────────────────────────────────────
def number_density(data_dir, tracer, z_eff):
    """n(z_eff) = (N_NGC + N_SGC) / (V_NGC + V_SGC), interpolated at z_eff."""
    data_dir = Path(data_dir)
    num = vol = 0.
    for cap in ('NGC', 'SGC'):
        gg = pd.read_csv(data_dir / 'nz' / f'{tracer}_{cap}_nz.txt', sep=r'\s+', header=None,
                         comment='#', names=['zmid', 'zlow', 'zhigh', 'nz', 'Nbin', 'Vol_bin'])
        num += interpolate.interp1d(gg['zmid'], gg['Nbin'], fill_value='extrapolate')(z_eff)
        vol += interpolate.interp1d(gg['zmid'], gg['Vol_bin'], fill_value='extrapolate')(z_eff)
    return float(num / vol)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Model errors (band around the model curve)
# ─────────────────────────────────────────────────────────────────────────────
# Number of small AbacusSummit boxes (500 h^-1 Mpc) behind each mock covariance
N_MOCKS = {'ELG1': 526, 'LRG3': 1299, 'ELGxLRG': 919}


def model_errors(data_dir, tracer, n_bins=26):
    """1-sigma model error: square root of the diagonal of the mock covariance,
    with the Hartlap factor (n - n_bins - 2) / (n - 1). Returns dict(e0, e2, ewp)."""
    cov_dir = Path(data_dir) / 'mock_covariances'
    n = N_MOCKS[tracer]
    hartlap = (n - n_bins - 2) / (n - 1)
    out = {}
    for key, stat in (('e0', 'mono'), ('e2', 'quad'), ('ewp', 'wp')):
        cov = np.loadtxt(cov_dir / f'{tracer}_cov_{stat}.txt', comments='#')
        out[key] = np.sqrt(np.diag(cov) * hartlap)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 4. chi^2 model vs DESI
# ─────────────────────────────────────────────────────────────────────────────
def chi2(model, data, data_dir, tracer, smin=None, smax=None):
    """chi^2 per statistic, for smin <= s (rp) <= smax:
         xi0, xi2 : full data covariance
         wp       : diagonal errors only (no wp covariance available)
    Returns dict(xi0, xi2, wp, nbins)."""
    data_dir = Path(data_dir)

    def mask(r):
        m = np.ones_like(r, dtype=bool)
        if smin is not None:
            m &= r >= smin
        if smax is not None:
            m &= r <= smax
        return m

    m = mask(data['s'])
    out = {}
    for key, stat in (('xi0', 'mono'), ('xi2', 'quad')):
        cov = np.loadtxt(data_dir / f'{tracer}_cov_{stat}.txt', comments='#')[np.ix_(m, m)]
        r = (data[key] - model[key])[m]
        out[key] = float(r @ np.linalg.solve(cov, r))        # r^T C^-1 r
    mp = mask(data['rp'])
    out['wp'] = float(np.sum(((data['wp'] - model['wp']) / data['ewp'])[mp]**2))
    out['nbins'] = int(m.sum())
    return out
