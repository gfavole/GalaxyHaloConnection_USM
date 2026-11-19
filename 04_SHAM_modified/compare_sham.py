# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Standard SHAM, modified SHAM (f_sat + V_peak threshold) and extended SHAM with the
Gaussian PDF of Favole+16, on the same (sub)halos and with the same seed, for
LRG3, ELG1 and the cross.

    python compare_sham.py      ->  outputs/figures/SHAM_compare_<T>.pdf
"""
from pathlib import Path

import config
from common import desi
from common.pipeline import run
from common.plots import plot_compare
from sham.mock import build_mock

RUNS = {'SHAM standard': ('SHAM_standard', config.SHAM_STANDARD),
        'SHAM f_sat + V_peak threshold': ('SHAM_modified', config.SHAM_MODIFIED),
        'SHAM Gaussian PDF (Favole+16)': ('SHAM_pdf', config.SHAM_PDF)}

if __name__ == '__main__':
    models = {}
    for label, (name, pars) in RUNS.items():
        config.SHAM = pars
        print(f"\n========== {label} ==========")
        models[label], _ = run(build_mock, name, config)

    figdir = Path(config.OUTDIR) / 'figures'
    for T in next(iter(models.values())):
        plot_compare({k: m[T] for k, m in models.items()}, desi.load_clustering(config.DATA_DIR, T), T,
                     figdir / f'SHAM_compare_{T}.pdf')
    print(f"Comparison figures: {figdir}/SHAM_compare_*.pdf")
