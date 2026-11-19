# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
The same HOD with two ways of placing the satellites:
  'particles'  on the DM particles of their host
  'nfw'        on an analytic NFW profile with Gaussian velocities
Same random seed, so centrals and number of satellites per halo are identical:
any difference in the clustering comes only from WHERE and with WHICH VELOCITY
the satellites are placed (one-halo term, Fingers of God in xi_2).

    python compare_placement.py      ->  outputs/figures/HOD_compare_placement_<T>.pdf
"""
from pathlib import Path

import config
from common import desi
from common.pipeline import run
from common.plots import plot_compare
from hod.mock import build_mock

if __name__ == '__main__':
    models = {}
    for placement in ('particles', 'nfw'):
        config.SAT_PLACEMENT = placement
        print(f"\n========== satellites: {placement} ==========")
        models[placement], _ = run(build_mock, f'HOD_{placement}', config)

    figdir = Path(config.OUTDIR) / 'figures'
    for T in models['nfw']:
        data = desi.load_clustering(config.DATA_DIR, T)
        plot_compare({f'HOD, sat. on {k}': m[T] for k, m in models.items()}, data, T,
                     figdir / f'HOD_compare_placement_{T}.pdf')
    print(f"Comparison figures: {figdir}/HOD_compare_placement_*.pdf")
