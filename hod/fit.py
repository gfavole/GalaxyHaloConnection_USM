# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
hod/fit.py
==========
Quick best fit of an HOD model to one DESI DR1 tracer (used by the fit_hod.py
scripts of 01_HOD_standard and 02_HOD_Rocher23).

At every step the halos are populated (same random seed every step, so the chi^2
surface is not noisy), xi0, xi2, wp are measured with FCFC and the chi^2 with
respect to DESI is computed (full covariance for xi0, xi2; diagonal for wp).
Minimisation with Nelder-Mead.

This is NOT a posterior: a few tens of FCFC runs that move the parameters to a
reasonable fit. Every step is appended to outputs/fit/fit_<model>_<tracer>.log;
at the end the best parameters are printed in config.py format and saved to
outputs/fit/fit_<model>_<tracer>.json.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from clustering.fcfc_utils import measure
from common import desi
from common.catalogs import read_halos, write_mock
from common.pipeline import n_target_of
from hod.mock import hosts_of, particles_of
from hod.model import populate, tune_density

# ── Free parameters, bounds and initial simplex step of each kind of fit ─────
ROCHER_FREE = dict(free=('logMc', 'sigma_M', 'logM0', 'alpha', 'f_sigma_v'), step=(0.15, 0.1, 0.3, 0.3, 0.2),
                   bounds=dict(logMc=(11.0, 13.0), sigma_M=(0.03, 1.0), logM0=(10.5, 13.5),
                               alpha=(-1.0, 2.0), f_sigma_v=(0.5, 3.0)))
FITS = {
    # Zheng: log M_min is set at every step by the density
    'zheng': dict(free=('sigma_logM', 'logM0', 'logM1', 'alpha'), step=(0.15, 0.4, 0.3, 0.2),
                  bounds=dict(sigma_logM=(0.05, 1.2), logM0=(10.5, 14.5), logM1=(12.0, 15.5),
                              alpha=(0.2, 2.0))),
    # Rocher+23 without conformity: A_c, A_s rescaled by a common factor (log M_1, gamma fixed)
    'rocher': ROCHER_FREE,
    # with conformity f_sat is fixed by the amplitudes, so it is a free parameter too
    'rocher_conf': dict(free=ROCHER_FREE['free'] + ('fsat',), step=ROCHER_FREE['step'] + (0.01,),
                        bounds={**ROCHER_FREE['bounds'], 'fsat': (0.002, 0.15)}),
}
LARGE = 5.0          # h^-1 Mpc: threshold of the "large scales" chi^2 printed at the end


def main(config, models, description=__doc__):
    """models: dict name -> (parameter dict of config, kind of fit in FITS, tracers allowed)."""
    ap = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tracer', choices=['LRG3', 'ELG1'])
    ap.add_argument('--model', choices=list(models), default=list(models)[0])
    ap.add_argument('--smin', type=float, default=getattr(config, 'SMIN', None),
                    help='fit only s, r_p >= smin (default: config SMIN, i.e. all scales)')
    ap.add_argument('--smax', type=float, default=getattr(config, 'SMAX', None))
    ap.add_argument('--stats', nargs='+', default=['xi0', 'xi2', 'wp'], choices=['xi0', 'xi2', 'wp'])
    ap.add_argument('--maxeval', type=int, default=30, help='maximum number of FCFC runs')
    args = ap.parse_args()
    T, M = args.tracer, args.model
    p0, kind, allowed = models[M]
    if T not in allowed:
        ap.error(f"model {M} is for {', '.join(allowed)}")
    p0 = dict(p0[T] if isinstance(p0, dict) and T in p0 else p0)
    FREE, BOUNDS, step = FITS[kind]['free'], FITS[kind]['bounds'], np.array(FITS[kind]['step'])
    tag = f'{M}_{T}'

    # ── set-up: data, halos, particles ──────────────────────────────────────
    outdir = Path(config.OUTDIR) / 'fit'
    outdir.mkdir(parents=True, exist_ok=True)
    log = open(outdir / f'fit_{tag}.log', 'a')
    L, fmt = config.BOXSIZE, getattr(config, 'CATALOG_FORMAT', 'hdf5')
    nthreads = getattr(config, 'FCFC_NTHREADS', 4)
    data_dir = Path(config.DATA_DIR)
    data = desi.load_clustering(data_dir, T)
    hosts = hosts_of(read_halos(config.HALO_CAT, config.COLUMNS, boxsize=L), config)
    _, n_target = n_target_of(config, T)
    placement = p0.get('placement', getattr(config, 'SAT_PLACEMENT', 'nfw'))
    particles = particles_of(config) if placement == 'particles' else None
    print(f"{T} ({M}): N_target = {n_target:,}, satellites on {placement}, fit {args.stats}, "
          f"s in [{args.smin or 'min'}, {args.smax or 'max'}] h^-1 Mpc, max {args.maxeval} FCFC runs")

    def evaluate(p):
        """Populate with parameters p, measure the clustering. Returns (p tuned, model, f_sat)."""
        p, _ = tune_density(hosts['logM'].values, p, n_target)
        rng = np.random.default_rng(getattr(config, 'SEED', 42))      # same seed every step
        mock = populate(hosts, p, config.Z_EFF, rng, placement=placement, particles=particles)
        f = write_mock(mock, outdir / f'fit_{tag}.{"h5" if fmt == "hdf5" else "txt"}', config.Z_EFF, config.OM, L, fmt)
        s, xi0, xi2, rp, wp = measure([f], outdir, config.FCFC_DIR, L, data['s'], fmt, nthreads, tag=f'fit_{tag}')
        fsat = float(mock['is_sat'].mean()) if 'is_sat' in mock else float('nan')
        return p, dict(s=s, xi0=xi0, xi2=xi2, rp=rp, wp=wp), fsat

    best = {'chi2': np.inf}
    neval = [0]

    def cost(x):
        """chi^2 of the parameter vector x (huge outside the bounds)."""
        p = {**p0, **dict(zip(FREE, x))}
        if any(not lo <= p[k] <= hi for k, (lo, hi) in BOUNDS.items()):
            return 1e30
        t = time.time()
        p, model, fsat = evaluate(p)
        c2 = desi.chi2(model, data, data_dir, T, args.smin, args.smax)
        tot = sum(c2[k] for k in args.stats)
        neval[0] += 1
        amp = f"logMmin {p['logMmin']:.3f}  " if kind == 'zheng' else f"Ac {p['Ac']:.4f} As {p['As']:.4f}  "
        line = (f"{neval[0]:3d}  chi2 = {tot:12.1f}  ( " + "  ".join(f"{k} {c2[k]:.1f}" for k in ('xi0', 'xi2', 'wp'))
                + " )  " + amp + "  ".join(f"{k} {p[k]:.3f}" for k in FREE)
                + f"  fsat {100 * fsat:.1f}%  [{time.time() - t:.0f} s]")
        print(line, flush=True)
        log.write(line + '\n'); log.flush()
        if tot < best['chi2']:
            best.update(chi2=tot, p=p, fsat=fsat)
        return tot

    # ── Nelder-Mead from the config values ──────────────────────────────────
    x0 = np.array([p0[k] for k in FREE])
    simplex = np.vstack([x0] + [x0 + np.eye(len(FREE))[i] * step[i] for i in range(len(FREE))])
    minimize(cost, x0, method='Nelder-Mead',
             options=dict(initial_simplex=simplex, maxfev=args.maxeval, xatol=0.02, fatol=1.0))

    # ── final report: chi^2 of the best model, all scales and large scales ──
    p, model, fsat = evaluate(best['p'])
    allc = desi.chi2(model, data, data_dir, T)
    large = desi.chi2(model, data, data_dir, T, LARGE, None)
    print(f"\nBest fit {T} ({M}) after {neval[0]} FCFC runs (f_sat = {100 * fsat:.1f}%):")
    print(f"{M.upper()} = dict(" + ", ".join(f"{k}={v!r}" if not isinstance(v, float) else f"{k}={v:.4g}"
                                         for k, v in p.items()) + ")")
    for name, c in (('all scales', allc), (f's >= {LARGE:g}', large)):
        print(f"  chi2 {name:12s}: xi0 {c['xi0']:.1f}  xi2 {c['xi2']:.1f}  wp {c['wp']:.1f}   [{c['nbins']} bins each]")
    out = dict(tracer=T, model=M, params=p, fsat=fsat, fit_stats=args.stats,
               smin=args.smin, smax=args.smax, n_eval=neval[0], chi2_all=allc, chi2_large=large)
    with open(outdir / f'fit_{tag}.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Saved {outdir / f'fit_{tag}.json'}")
