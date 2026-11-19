# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Standard HOD (Zheng et al. 2007) on the AbacusSummit host halos, compared with
DESI DR1 LRG3, ELG1 and their cross-correlation.

    python run_hod.py

Outputs in outputs/: mocks HOD_<placement>_<T>.h5, model clustering
HOD_<placement>_<T>_clustering.txt, figures in outputs/figures/ and a summary json
(number of galaxies, satellite fraction, chi^2).
"""
import config                                   # this folder's config.py (also sets the paths)
from common.pipeline import run
from hod.mock import build_mock

if __name__ == '__main__':
    run(build_mock, f'HOD_{config.SAT_PLACEMENT}', config)
