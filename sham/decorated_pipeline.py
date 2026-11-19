# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
sham/decorated_pipeline.py
==========================
Workflow of the decorated SHAM on IllustrisTNG (05_SHAM_decorated). Unlike the other
models, the comparison is not with DESI but with the "truth" of the hydrodynamical
simulation, galaxy by galaxy selected in the same way:

  1. read the hydro galaxies (TNG300-1) and the dark subhalos (TNG300-1-Dark)
  2. truth: LRG-like and ELG-like samples of the hydro galaxies at the DESI DR1 densities
  3. model: standard SHAM on V_peak of the dark subhalos + decoration with a secondary
     subhalo property (sham/decorated.py), same selections
  4. FCFC: xi0, xi2, wp of truth and models (and the ELG x LRG cross)
  5. figures: models vs truth (outputs/figures/)
"""

import json
from pathlib import Path

import numpy as np

from clustering.fcfc_utils import measure
from common import desi
from common.catalogs import read_halos, write_mock
from common.plots import plot_compare
from sham.decorated import decorate, select_tracer, sham_stellar_mass

GAL_COLUMNS = dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz", is_sat="is_sat",
                   logMstar="logMstar", log_sSFR="log_sSFR", g_i="g_i")
SUB_COLUMNS = dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz", is_sat="is_sat",
                   Vpeak="Vpeak", z_starve="z_starve", c_infall="c_infall")


def simulation_attrs(path):
    """Box size, redshift and Omega_m stored as attributes of the TNG catalog."""
    import h5py
    with h5py.File(path, 'r') as f:
        return float(f.attrs['boxsize']), float(f.attrs['redshift']), float(f.attrs['Om'])


def load(config):
    """Hydro galaxies, dark subhalos, box size, redshift, Omega_m."""
    L, z, om = simulation_attrs(config.TNG_HYDRO_CAT)
    hydro = read_halos(config.TNG_HYDRO_CAT, GAL_COLUMNS, boxsize=L)
    dark = read_halos(config.TNG_DARK_CAT, SUB_COLUMNS, boxsize=L)
    dark = dark[np.isfinite(dark['Vpeak'].values)].reset_index(drop=True)
    print(f"    hydro: {len(hydro):,} galaxies;  dark: {len(dark):,} subhalos;  box {L} h^-1 Mpc, z = {z:.2f}")
    return hydro, dark, L, z, om


def targets(config, L):
    """Number of galaxies of each sample: DESI DR1 density at z_eff times the box volume."""
    return {T: int(round(desi.number_density(config.DATA_DIR, T, config.Z_EFF) * L**3))
            for T in config.TRACERS}


def build_model(dark, hydro, p, n_target, rng):
    """Decorated SHAM mock for every tracer: dict tracer -> DataFrame."""
    mock = dark.copy()
    mock['logMstar'] = sham_stellar_mass(mock['Vpeak'].values, hydro['logMstar'].values,
                                         p['sigma_logMstar'], rng)
    mock = mock[np.isfinite(mock['logMstar'].values)].reset_index(drop=True)
    mock = decorate(mock, hydro, p['secondary'], p['rho'], p['mstar_bins'], rng)
    return {T: select_tracer(mock, T, n) for T, n in n_target.items()}


def clustering(samples, name, outdir, config, L, z, om, scale):
    """FCFC for each sample and the cross. samples: dict tracer -> DataFrame."""
    fmt, nth = getattr(config, 'CATALOG_FORMAT', 'hdf5'), getattr(config, 'FCFC_NTHREADS', 4)
    files = {T: write_mock(s, outdir / f'{name}_{T}.{"h5" if fmt == "hdf5" else "txt"}', z, om, L, fmt)
             for T, s in samples.items()}
    out = {}
    for T, f in files.items():
        s, xi0, xi2, rp, wp = measure([f], outdir, config.FCFC_DIR, L, scale, fmt, nth, tag=f'{name}_{T}')
        out[T] = dict(s=s, xi0=xi0, xi2=xi2, rp=rp, wp=wp)
    if {'ELG1', 'LRG3'} <= set(files):
        s, xi0, xi2, rp, wp = measure([files['ELG1'], files['LRG3']], outdir, config.FCFC_DIR, L, scale, fmt,
                                      nth, tag=f'{name}_ELGxLRG')
        out['ELGxLRG'] = dict(s=s, xi0=xi0, xi2=xi2, rp=rp, wp=wp)
    return out


def summary(samples):
    """Satellite fraction and median M_* of each sample."""
    return {T: dict(N=len(s), fsat=float(s['is_sat'].mean()), median_logMstar=float(np.median(s['logMstar'])))
            for T, s in samples.items()}


def run(config, models, fig_tag='SHAM_decorated'):
    """models: dict label -> parameter dict. Runs truth + models and the comparison figures
    outputs/figures/<fig_tag>_vs_TNG_<T>.pdf."""
    outdir = Path(config.OUTDIR)
    (outdir / 'figures').mkdir(parents=True, exist_ok=True)
    scale = desi.load_clustering(config.DATA_DIR, 'LRG3')['s']      # the same 26 bins as the course

    print("[1] TNG catalogs")
    hydro, dark, L, z, om = load(config)
    n_target = targets(config, L)
    print(f"    N per sample: {n_target}")

    print("[2] truth: hydrodynamical galaxies")
    truth_samples = {T: select_tracer(hydro, T, n) for T, n in n_target.items()}
    report = {'truth': summary(truth_samples)}
    truth = clustering(truth_samples, 'TNG_hydro', outdir, config, L, z, om, scale)

    print("[3] decorated SHAM on the dark-matter-only simulation")
    results = {}
    for label, p in models.items():
        rng = np.random.default_rng(getattr(config, 'SEED', 42))
        samples = build_model(dark, hydro, p, n_target, rng)
        report[label] = summary(samples)
        tag = 'SHAM_' + ''.join(c if c.isalnum() else '_' for c in label)
        results[label] = clustering(samples, tag, outdir, config, L, z, om, scale)
        print(f"    {label}: " + ", ".join(f"{T} f_sat = {100 * r['fsat']:.1f}%" for T, r in report[label].items()))
    print("    truth: " + ", ".join(f"{T} f_sat = {100 * r['fsat']:.1f}%" for T, r in report['truth'].items()))

    print("[4] figures")
    for T, tr in truth.items():
        data = dict(tr, e0=np.zeros_like(tr['xi0']), e2=np.zeros_like(tr['xi2']), ewp=np.zeros_like(tr['wp']))
        path = outdir / 'figures' / f'{fig_tag}_vs_TNG_{T}.pdf'
        plot_compare({k: r[T] for k, r in results.items()}, data, T, path, data_label=f'TNG300-1 hydro {T}')
        # large-scale ratio model / truth (mean over s > 5 h^-1 Mpc)
        big = tr['s'] > 5.
        for k, r in results.items():
            report[k][T + '_xi0_ratio_large'] = float(np.mean(r[T]['xi0'][big] / tr['xi0'][big]))
    with open(outdir / f'{fig_tag}_summary.json', 'w') as f:
        json.dump(report, f, indent=2)
    print(f"Figures: {outdir}/figures/{fig_tag}_vs_TNG_*.pdf;  summary: {fig_tag}_summary.json")
    return truth, results, report
