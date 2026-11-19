# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Modified HOD of Rocher et al. 2023 for the ELGs (config.ELG_MODEL), with the LRGs
on the standard Zheng HOD, compared with DESI DR1 LRG3, ELG1 and ELG x LRG.

    python run_rocher.py
"""
import config
from common.pipeline import run
from hod.mock import build_mock

if __name__ == '__main__':
    p = config.ELG_MODEL
    name = (f"HOD_Rocher23_{p['model']}" + ('_vbias' if p.get('f_sigma_v', 1.) != 1. else '')
            + ('_conf' if p.get('conformity') else ''))
    run(build_mock, name, config)
