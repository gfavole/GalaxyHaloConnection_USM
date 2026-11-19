# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
One modified SHAM (config.SHAM: SHAM_MODIFIED or SHAM_PDF) compared with DESI DR1
LRG3, ELG1 and ELG x LRG.

    python run_sham.py
"""
import config
from common.pipeline import run
from sham.mock import build_mock

if __name__ == '__main__':
    model = config.SHAM['ELG1'].get('model', 'standard')
    run(build_mock, f'SHAM_{model}', config)
