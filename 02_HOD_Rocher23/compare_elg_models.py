# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
ELG1 with the HOD models of config.ELG_MODELS (same halos, same seed):
  Zheng (complete sample) -> Gaussian centrals + velocity bias
  -> + strict conformity + modified satellite profile
plus, if CROSS is True, the ELG x LRG cross with the LRG Zheng model.

    python compare_elg_models.py   ->  outputs/figures/HOD_compare_ELGmodels_{ELG1,ELGxLRG}.pdf
"""
import re
from pathlib import Path

import config
from common import desi
from common.pipeline import run
from common.plots import plot_compare
from hod.mock import build_mock

if __name__ == '__main__':
    models, fixed = {}, None
    for label, p in config.ELG_MODELS.items():
        config.HOD = {**config.HOD, 'ELG1': p}                  # only the ELG model changes
        tag = 'HOD_ELG_' + re.sub(r'[^A-Za-z0-9]+', '_', label).strip('_')
        print(f"\n========== ELG1: {label} ==========")
        models[label], summary = run(build_mock, tag, config, fixed=fixed)
        if fixed is None and 'LRG3' in summary:
            # the LRGs (Zheng) are the same in every run: build and measure them only once
            fixed = {'LRG3': {**summary['LRG3'], 'model': models[label]['LRG3']}}

    figdir = Path(config.OUTDIR) / 'figures'
    for T in ('ELG1', 'ELGxLRG'):
        if all(T in m for m in models.values()):
            plot_compare({k: m[T] for k, m in models.items()}, desi.load_clustering(config.DATA_DIR, T), T,
                         figdir / f'HOD_compare_ELGmodels_{T}.pdf')
    print(f"Comparison figures: {figdir}/HOD_compare_ELGmodels_*.pdf")
