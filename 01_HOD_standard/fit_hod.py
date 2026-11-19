# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Quick best fit of the standard HOD (Zheng+07) to one DESI DR1 tracer.
Free: sigma_logM, log M_0, log M_1, alpha; log M_min is set at every step by the
DESI number density. Nelder-Mead, a few tens of FCFC runs (see hod/fit.py).

    python fit_hod.py LRG3                          # all scales, xi0 + xi2 + wp
    python fit_hod.py ELG1 --smin 5                 # only s, r_p >= 5 h^-1 Mpc
    python fit_hod.py LRG3 --stats xi0 wp           # chi^2 of xi0 and wp only

Paste the printed best fit into config.py and run run_hod.py for the figures.
"""
import config
from hod.fit import main

# name -> (parameters in config, kind of fit, tracers)
MODELS = {'zheng': (config.HOD, 'zheng', ('LRG3', 'ELG1'))}

if __name__ == '__main__':
    main(config, MODELS, description=__doc__)
