# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Standard SHAM on V_peak (hosts -> centrals, satellite candidates -> satellites),
compared with DESI DR1 LRG3, ELG1 and ELG x LRG.

    python run_sham.py

Outputs in outputs/: mocks SHAM_standard_<T>.h5 (with log M_* and V_peak),
clustering SHAM_standard_<T>_clustering.txt, figures and a summary json.
"""
import config
from common.pipeline import run
from sham.mock import build_mock

if __name__ == '__main__':
    run(build_mock, 'SHAM_standard', config)
