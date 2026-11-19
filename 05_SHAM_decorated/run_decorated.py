# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Decorated SHAM on TNG300-1-Dark vs the true galaxies of TNG300-1 (hydro):
the models of config.MODELS (SHAM without secondary property, with z_starve and
with c_infall), LRG-like and ELG-like samples at the DESI DR1 densities.

    python run_decorated.py

Outputs in outputs/: mocks and truth samples (.h5), figures
figures/SHAM_decorated_vs_TNG_{LRG3,ELG1,ELGxLRG}.pdf and SHAM_decorated_summary.json
(satellite fractions, median M_*, large-scale ratio model / truth).
"""
import config
from sham.decorated_pipeline import run

if __name__ == '__main__':
    run(config, config.MODELS)
