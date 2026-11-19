# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Quick best fit of the Rocher+23 ELG models to DESI DR1 ELG1 (see hod/fit.py).

  --model mhmq_vb        asymmetric Gaussian (mHMQ) + velocity bias (config.MHMQ_VB)
  --model mhmq_vb_conf   + strict conformity + modified profile (config.MHMQ_VB_CONF)
  --model ghod_vb        symmetric Gaussian (GHOD, base model) + velocity bias (config.GHOD_VB)
  --model zheng          Zheng for LRG3 or ELG1 (as in 01_HOD_standard)

Free: log M_c, sigma_M, log M_0, alpha, f_sigma_v (and f_sat with conformity);
A_c, A_s are set by the DESI density (and f_sat); log M_1 and gamma are fixed.

    python fit_hod.py ELG1 --model ghod_vb
    python fit_hod.py ELG1 --model mhmq_vb
    python fit_hod.py ELG1 --model mhmq_vb_conf --maxeval 20
"""
import config
from hod.fit import main

# name -> (parameters in config, kind of fit, tracers)
MODELS = {'mhmq_vb': (config.MHMQ_VB, 'rocher', ('ELG1',)),
          'mhmq_vb_conf': (config.MHMQ_VB_CONF, 'rocher_conf', ('ELG1',)),
          'ghod_vb': (config.GHOD_VB, 'rocher', ('ELG1',)),
          'zheng': ({'LRG3': config.ZHENG_LRG, 'ELG1': config.ZHENG_ELG}, 'zheng', ('LRG3', 'ELG1'))}

if __name__ == '__main__':
    main(config, MODELS, description=__doc__)
