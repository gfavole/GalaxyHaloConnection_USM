# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
common/plots.py
===============
Model vs DESI: s xi0, s xi2 and rp wp in a three-panel figure.
Colours: LRG red, ELG blue, ELG x LRG cross black; DESI data as black points.
"""

import matplotlib
matplotlib.use('Agg')                    # no window: save straight to PDF
import matplotlib.pyplot as plt

COLOR = {'LRG3': '#B83227', 'ELG1': '#2F6DB5', 'ELGxLRG': 'k'}

# (x axis, statistic, error, y label, x label)
PANELS = [('s', 'xi0', 'e0', r'$s\,\xi_0(s)$', r'$s\ [h^{-1}{\rm Mpc}]$'),
          ('s', 'xi2', 'e2', r'$s\,\xi_2(s)$', r'$s\ [h^{-1}{\rm Mpc}]$'),
          ('rp', 'wp', 'ewp', r'$r_p\,w_p(r_p)$', r'$r_p\ [h^{-1}{\rm Mpc}]$')]


def _data(a, data, r, y, e, tracer, label=None):
    """Data points with error bars (DESI by default)."""
    if data is not None:
        a.errorbar(data[r], data[r] * data[y], data[r] * data[e], color='k', ls='none',
                   marker='o', ms=4, capsize=3, label=label or f'DESI DR1 {tracer}')


def plot_clustering(model, data, tracer, model_label, outpath, model_err=None):
    """One model. model_err: dict(e0, e2, ewp) from the mock covariances (shaded band)."""
    c = COLOR.get(tracer, 'k')
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for a, (r, y, e, ylab, xlab) in zip(ax, PANELS):
        _data(a, data, r, y, e, tracer)
        a.plot(model[r], model[r] * model[y], color=c, lw=2, label=model_label)
        if model_err is not None:
            a.fill_between(model[r], model[r] * (model[y] - model_err[e]),
                           model[r] * (model[y] + model_err[e]), color=c, alpha=0.2, lw=0)
        a.set_xscale('log'); a.set_xlabel(xlab); a.set_ylabel(ylab)
    ax[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    return outpath


def plot_compare(models, data, tracer, outpath, data_label=None):
    """Several models of the same tracer: same colour, different line styles."""
    c = COLOR.get(tracer, 'k')
    styles = ['-', '--', ':', '-.', (0, (5, 1, 1, 1))]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for a, (r, y, e, ylab, xlab) in zip(ax, PANELS):
        _data(a, data, r, y, e, tracer, data_label)
        for (lab, m), ls in zip(models.items(), styles):
            a.plot(m[r], m[r] * m[y], color=c, lw=2, ls=ls, label=lab)
        a.set_xscale('log'); a.set_xlabel(xlab); a.set_ylabel(ylab)
    ax[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    return outpath
