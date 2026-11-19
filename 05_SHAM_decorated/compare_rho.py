# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
How much assembly bias do we need? The decorated SHAM with z_starve and different
correlation strengths rho (0 = no assembly bias, 1 = full rank order), compared with
TNG300-1.

    python compare_rho.py       ->  outputs/figures/SHAM_decorated_rho_vs_TNG_<T>.pdf
"""
import config
from sham.decorated_pipeline import run

RHO_VALUES = (0.0, 0.5, 1.0)

if __name__ == '__main__':
    run(config, {f"z_starve, rho = {r}": {**config.DECORATED, "rho": r} for r in RHO_VALUES},
        fig_tag='SHAM_decorated_rho')
