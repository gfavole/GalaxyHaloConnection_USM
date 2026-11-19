# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
sham/mock.py
============
build_mock(halos, tracer, n_target, config, rng) for common/pipeline.py.

The parameters of each tracer are in config.SHAM[tracer]; the key 'model' chooses
how the galaxy sample is selected after the abundance matching:

  'standard'   the n_target objects with the largest V'_peak (or M_* >= logMstar_min)
  'modified'   fixed satellite fraction 'fsat' and/or random draw above 'Vpeak_min'
               (satellites above 'Vpeak_min_sat' if given)
  'pdf'        Gaussian PDF in V'_peak (V_mean, sigma_V) with fixed 'fsat' (Favole+16)

The course catalogs come in two files: the host halos (config.HALO_CAT, centrals)
and the satellite candidates (config.PARTICLE_CAT: DM particles with a sampled
V_peak). Both are matched together: hosts get is_sat = 0, candidates is_sat = 1.
A single catalog with its own is_sat column (PARTICLE_CAT = None) also works.
"""

import pandas as pd

from common.catalogs import read_halos
from common.pipeline import host_fraction_of
from sham.model import abundance_matching, select
from sham.pdf import select_pdf

_OBJECTS = {}          # cache: hosts + candidates are read only once per run


def objects_of(halos, config):
    """Hosts + satellite candidates in one table, read once."""
    if 'df' not in _OBJECTS:
        part_cat = getattr(config, 'PARTICLE_CAT', None)
        if part_cat:
            cols = getattr(config, 'SHAM_PARTICLE_COLUMNS',
                           dict(x="x", y="y", z="z", vx="vx", vy="vy", vz="vz", Vpeak="Vpeak"))
            sats = read_halos(part_cat, cols, boxsize=config.BOXSIZE)
            print(f"    {len(sats):,} satellite candidates read")
            hosts = halos.copy()
            hosts['is_sat'] = 0
            sats['is_sat'] = 1
            keep = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'Vpeak', 'is_sat']
            extra = [c for c in halos.columns if c not in keep + ['host_id']]   # e.g. conc, logM
            _OBJECTS['df'] = pd.concat([hosts[keep + extra], sats[keep]], ignore_index=True)
        else:
            _OBJECTS['df'] = halos
    return _OBJECTS['df']


def build_mock(halos, tracer, n_target, config, rng):
    p = dict(config.SHAM[tracer])
    model = p.get('model', 'standard')
    objs = objects_of(halos, config)

    # 1. abundance matching. The catalogs keep a random fraction f of the halos:
    #    matching in the effective volume f L^3 keeps the M_* of each V_peak as in the full box.
    L_eff = config.BOXSIZE * host_fraction_of(config) ** (1. / 3.)
    am = abundance_matching(objs, p.get('sigma_AM', 0.), L_eff, rng)

    # 2. selection of the sample
    if model == 'standard':
        if p.get('logMstar_min') is not None:
            return select(am, logMstar_min=p['logMstar_min'])
        return select(am, n_target=n_target)
    if model == 'modified':
        return select(am, n_target=n_target, fsat=p.get('fsat'), vpeak_min=p.get('Vpeak_min'), rng=rng,
                      vpeak_min_sat=p.get('Vpeak_min_sat'))
    if model == 'pdf':
        return select_pdf(am, n_target, p['fsat'], p['V_mean'], p['sigma_V'], rng,
                          v_mean_sat=p.get('V_mean_sat'), sigma_v_sat=p.get('sigma_V_sat'))
    raise ValueError("SHAM model must be 'standard', 'modified' or 'pdf' (decorated: sham/decorated_pipeline.py)")
