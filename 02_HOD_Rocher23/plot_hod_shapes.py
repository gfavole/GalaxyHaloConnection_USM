# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Mean halo occupation <N_cen>(M), <N_sat>(M) of the ELG models in config.ELG_MODELS,
after matching the DESI density on the host-halo catalog (same tuning as the runs).
With conformity the satellite curve is the mean over all halos, <N_cen> <N_sat|cen>.
No FCFC needed: it only reads the halo masses.

    python plot_hod_shapes.py      ->  outputs/figures/HOD_shapes_ELG1.pdf
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import config
from common.catalogs import read_halos
from common.pipeline import n_target_of
from hod.mock import hosts_of
from hod.model import mean_occupations, tune_density

# same line styles, in the same order, as the clustering comparison (common/plots.plot_compare)
STYLES = [dict(ls=ls, lw=2) for ls in ('-', '--', ':', '-.', (0, (5, 1, 1, 1)))]


def occupations(p, logM):
    """<N_cen>, and <N_sat> averaged over all halos (also with conformity)."""
    nc, ns = mean_occupations(logM, p)
    return nc, (nc * ns if p.get('conformity', False) else ns)


def plot_shapes(models, outpath, tracer='ELG1', color='#2F6DB5'):
    """models: dict label -> parameter dict already tuned to the density."""
    logM = np.linspace(10.5, 15.0, 400)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    for (label, p), st in zip(models.items(), STYLES):
        nc, ns = occupations(p, logM)
        axes[0].plot(logM, nc, color=color, label=label, **st)
        axes[1].plot(logM, ns, color=color, label=label, **st)
    for a, t in zip(axes, (r'$\langle N_{\rm cen}\rangle$', r'$\langle N_{\rm sat}\rangle$')):
        a.set_yscale('log'); a.set_ylim(1e-4, 3); a.set_xlim(10.5, 15)
        a.set_xlabel(r'$\log_{10} M_{\rm h}\ [h^{-1}M_\odot]$'); a.set_title(t)
    axes[0].set_ylabel('mean occupation')
    # one legend under both panels, so it does not cover the curves
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=10)
    fig.suptitle(f'{tracer}: shape of the HOD', y=0.99)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(outpath)
    print(f"Saved {outpath}")


if __name__ == '__main__':
    figdir = Path(config.OUTDIR) / 'figures'
    figdir.mkdir(parents=True, exist_ok=True)
    hosts = hosts_of(read_halos(config.HALO_CAT, config.COLUMNS, boxsize=config.BOXSIZE), config)
    _, n_target = n_target_of(config, 'ELG1')
    tuned = {}
    for label, p in config.ELG_MODELS.items():
        tuned[label], msg = tune_density(hosts['logM'].values, dict(p), n_target)
        print(f"{label}: {msg}")
    plot_shapes(tuned, figdir / 'HOD_shapes_ELG1.pdf')
