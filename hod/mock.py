# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
hod/mock.py
===========
build_mock(halos, tracer, n_target, config, rng) for common/pipeline.py.

The HOD parameters of each tracer are in config.HOD[tracer]. If config.TUNE_DENSITY
is True, the expected number of galaxies is matched to the DESI number density
(times the box volume) by adjusting log M_min (Zheng) or the amplitudes A_c, A_s
(GHOD, mHMQ). Satellites go on the DM particles of the host or on an NFW profile
(config.SAT_PLACEMENT, or the 'placement' key of the model).
"""

from common.catalogs import read_halos
from hod.model import populate, tune_density

_PARTICLES = {}            # cache: the particle catalog is read only once per run


def particles_of(config):
    """DM particles (satellite candidates), read once."""
    if 'df' not in _PARTICLES:
        df = read_halos(config.PARTICLE_CAT, config.PARTICLE_COLUMNS, boxsize=config.BOXSIZE)
        if 'vxhost' in df:                       # catalogs with velocities relative to the host
            for c in ('vx', 'vy', 'vz'):
                df[c] = df[c].values + df[c + 'host'].values
            df = df.drop(columns=['vxhost', 'vyhost', 'vzhost'])
        _PARTICLES['df'] = df
        print(f"    {len(_PARTICLES['df']):,} DM particles read")
    return _PARTICLES['df']


def hosts_of(halos, config):
    """Host halos with Rvir in h^-1 Mpc."""
    hosts = halos[halos['is_sat'] == 0] if 'is_sat' in halos else halos
    if 'Rvir' in hosts and getattr(config, 'RVIR_TO_MPC', 1.) != 1.:
        hosts = hosts.assign(Rvir=hosts['Rvir'].values * config.RVIR_TO_MPC)
    return hosts


def build_mock(halos, tracer, n_target, config, rng):
    p = dict(config.HOD[tracer])
    hosts = hosts_of(halos, config)
    # 1. density: adjust log M_min or the amplitudes to get n_target galaxies
    if getattr(config, 'TUNE_DENSITY', True):
        p, msg = tune_density(hosts['logM'].values, p, n_target)
        print(f"    {tracer} ({p.get('model', 'zheng')}{', conformity' if p.get('conformity') else ''}): {msg}")
    # 2. where to put the satellites
    placement = p.get('placement', getattr(config, 'SAT_PLACEMENT', 'nfw'))
    particles = particles_of(config) if placement == 'particles' else None
    # 3. populate
    return populate(hosts, p, config.Z_EFF, rng, placement=placement, particles=particles)
