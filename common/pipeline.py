# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
common/pipeline.py
==================
The same workflow for every model of the course:

  1. DESI DR1 data (xi0, xi2, wp), covariances and n(z_eff)
  2. host-halo catalog
  3. one mock per tracer (LRG3, ELG1), built by the model, with the DESI number density
  4. FCFC: xi0, xi2, wp of each mock and, if CROSS = True, the ELG x LRG cross-correlation
  5. chi^2 with respect to DESI, figures (outputs/figures/) and text files with the model (outputs/)

A model only has to provide one function

    build_mock(halos, tracer, n_target, config, rng) -> DataFrame

with columns x, y, z, vx, vy, vz (and optionally is_sat, logM, logMstar, ...).
Changing model means changing one function: everything else is identical.
"""

import json
import time
from pathlib import Path

import numpy as np

from clustering.fcfc_utils import measure
from common import desi
from common.catalogs import read_halos, write_mock
from common.plots import plot_clustering


def host_fraction_of(config):
    """Random fraction of halos kept in the host catalog (to make it lighter): the
    'host_fraction' attribute of the HDF5 file if present, else config.HOST_FRACTION.
    The mock then has N = n_DESI L^3 host_fraction galaxies: the same galaxies as with
    all the halos, randomly diluted, so the clustering does not change."""
    try:
        import h5py
        with h5py.File(config.HALO_CAT, 'r') as f:
            if 'host_fraction' in f.attrs:
                return float(f.attrs['host_fraction'])
    except (OSError, ImportError):
        pass
    return float(getattr(config, 'HOST_FRACTION', 1.0))


def n_target_of(config, tracer):
    """Number of galaxies of the mock: n_DESI(z_eff) L^3 host_fraction."""
    n = desi.number_density(config.DATA_DIR, tracer, config.Z_EFF)
    return n, int(round(n * config.BOXSIZE**3 * host_fraction_of(config)))


_HALOS = {}                                       # halo catalogs already read, by path


def run(build_mock, model_name, config, fixed=None):
    """Run the full model. Returns (model measurements, summary).

    fixed : optional {tracer: summary entry of an earlier run} for tracers that do
            not change between runs (e.g. the LRGs when only the ELG model changes):
            their mock and clustering are reused instead of being rebuilt."""
    fixed = fixed or {}
    t0 = time.time()
    data_dir = Path(config.DATA_DIR)
    outdir = Path(config.OUTDIR)
    figdir = outdir / 'figures'
    figdir.mkdir(parents=True, exist_ok=True)
    L = config.BOXSIZE
    fmt = getattr(config, 'CATALOG_FORMAT', 'hdf5')
    nthreads = getattr(config, 'FCFC_NTHREADS', 4)
    rng = np.random.default_rng(getattr(config, 'SEED', 42))   # fixed seed: reproducible results
    tracers = list(getattr(config, 'TRACERS', ('LRG3', 'ELG1')))

    # ── 1. Data ─────────────────────────────────────────────────────────────
    print(f"\n[1] DESI DR1 data ({', '.join(tracers)})")
    data = {T: desi.load_clustering(data_dir, T) for T in tracers}
    if getattr(config, 'CROSS', False):
        data['ELGxLRG'] = desi.load_clustering(data_dir, 'ELGxLRG')
    scale = data[tracers[0]]['s']                # the model is interpolated at the s of the data

    # ── 2. Halos ────────────────────────────────────────────────────────────
    print("[2] Halo catalog")
    if config.HALO_CAT not in _HALOS:              # read once per session (83M rows for Abacus)
        _HALOS.clear()
        _HALOS[config.HALO_CAT] = read_halos(config.HALO_CAT, config.COLUMNS, boxsize=L)
    halos = _HALOS[config.HALO_CAT]
    print(f"    {len(halos):,} halos, box {L} h^-1 Mpc, host fraction {host_fraction_of(config)}")

    # ── 3. Mocks ────────────────────────────────────────────────────────────
    print("[3] Mocks")
    mock_files, summary = {}, {}
    for T in tracers:
        if T in fixed:                             # unchanged tracer: reuse the earlier mock
            mock_files[T] = Path(fixed[T]['mock_file'])
            summary[T] = {k: v for k, v in fixed[T].items() if k not in ('chi2', 'model')}
            print(f"    {T}: reused from the previous run (N_mock = {summary[T]['N_mock']:,})")
            continue
        n, n_target = n_target_of(config, T)
        mock = build_mock(halos, T, n_target, config, rng)      # <-- the model enters here
        fsat = float(mock['is_sat'].mean()) if 'is_sat' in mock else float('nan')
        print(f"    {T}: n(z_eff) = {n:.3e} h^3 Mpc^-3 -> N_target = {n_target:,};"
              f" N_mock = {len(mock):,}; f_sat = {100 * fsat:.1f}%")
        mock_files[T] = write_mock(mock, outdir / f'{model_name}_{T}.{"h5" if fmt == "hdf5" else "txt"}',
                                   config.Z_EFF, config.OM, L, fmt)
        summary[T] = dict(n_data=n, N_target=n_target, N_mock=len(mock), fsat=fsat,
                          mock_file=str(mock_files[T]))

    # ── 4. Clustering with FCFC ─────────────────────────────────────────────
    print("[4] FCFC")
    model = {}
    for T in tracers:
        if T in fixed and 'model' in fixed[T]:    # unchanged tracer: reuse its clustering
            model[T] = fixed[T]['model']
            continue
        s, xi0, xi2, rp, wp = measure([mock_files[T]], outdir, config.FCFC_DIR, L, scale, fmt, nthreads,
                                      tag=f'{model_name}_{T}')
        model[T] = dict(s=s, xi0=xi0, xi2=xi2, rp=rp, wp=wp)
    if getattr(config, 'CROSS', False) and {'ELG1', 'LRG3'} <= set(tracers):
        s, xi0, xi2, rp, wp = measure([mock_files['ELG1'], mock_files['LRG3']], outdir, config.FCFC_DIR, L,
                                      scale, fmt, nthreads, tag=f'{model_name}_ELGxLRG')
        model['ELGxLRG'] = dict(s=s, xi0=xi0, xi2=xi2, rp=rp, wp=wp)

    # ── 5. chi^2 and figures ────────────────────────────────────────────────
    print("[5] chi^2 and figures")
    for T, m in model.items():
        np.savetxt(outdir / f'{model_name}_{T}_clustering.txt',
                   np.column_stack([m['s'], m['xi0'], m['xi2'], m['wp']]),
                   header='s(=rp) xi0 xi2 wp   [h^-1 Mpc]')
        c2 = desi.chi2(m, data[T], data_dir, T, getattr(config, 'SMIN', None), getattr(config, 'SMAX', None))
        summary.setdefault(T, {})['chi2'] = c2
        print(f"    {T}: chi2(xi0) = {c2['xi0']:.1f}  chi2(xi2) = {c2['xi2']:.1f}  "
              f"chi2(wp) = {c2['wp']:.1f}   [{c2['nbins']} bins each]")
        plot_clustering(m, data[T], T, model_name, figdir / f'{model_name}_{T}.pdf',
                        model_err=desi.model_errors(data_dir, T))

    with open(outdir / f'{model_name}_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nDone in {time.time() - t0:.0f} s. Outputs in {outdir}")
    return model, summary
