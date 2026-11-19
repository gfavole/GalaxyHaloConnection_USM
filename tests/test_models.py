# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Quick tests of the models, without FCFC and without the course catalogs:

    python -m pytest -q tests
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.desi import number_density                                  # noqa: E402
from hod.model import (expected_number, ncen_mean, nsat_mean, populate,  # noqa: E402
                       sample_nfw_radius, tune_density, tune_logMmin)
from sham.decorated import decorate, select_tracer, sham_stellar_mass   # noqa: E402
from sham.model import abundance_matching, select, weighted_draw         # noqa: E402
from sham.pdf import select_pdf                                          # noqa: E402

DATA = Path(__file__).resolve().parents[1] / 'data' / 'desi_dr1'
ZHENG = dict(model='zheng', logMmin=12.5, sigma_logM=0.3, logM0=12.5, logM1=13.5, alpha=1.0)
GHOD = dict(model='ghod', Ac=0.1, logMc=11.9, sigma_M=0.2, As=0.01, logM0=11.7, logM1=13., alpha=0.5,
            placement='nfw')


def halos(n=200000, seed=0):
    """Fake host halos uniform in a 500 h^-1 Mpc box."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame(dict(x=rng.uniform(0, 500, n), y=rng.uniform(0, 500, n), z=rng.uniform(0, 500, n),
                             vx=np.zeros(n), vy=np.zeros(n), vz=np.zeros(n),
                             logM=rng.uniform(11, 15, n), Rvir=np.full(n, 0.3)))


def objects(n_host=100000, n_sat=300000, seed=0):
    """Fake hosts + satellite candidates with V_peak and a concentration."""
    rng = np.random.default_rng(seed)
    v = np.r_[rng.lognormal(np.log(200), 0.4, n_host), rng.lognormal(np.log(80), 0.4, n_sat)]
    n = n_host + n_sat
    return pd.DataFrame(dict(x=rng.uniform(0, 500, n), y=rng.uniform(0, 500, n), z=rng.uniform(0, 500, n),
                             vx=np.zeros(n), vy=np.zeros(n), vz=np.zeros(n), Vpeak=v,
                             conc=np.r_[rng.lognormal(np.log(8), 0.3, n_host), np.full(n_sat, np.nan)],
                             is_sat=np.r_[np.zeros(n_host, int), np.ones(n_sat, int)]))


# ── HOD ─────────────────────────────────────────────────────────────────────
def test_zheng_occupations():
    assert np.isclose(ncen_mean(12.5, 12.5, 0.3), 0.5)            # <N_cen> = 1/2 at M_min
    assert nsat_mean(12.0, 12.5, 0.3, 12.5, 13.5, 1.0) == 0.      # no satellites below M_0


def test_populate_matches_expectation():
    h = halos()
    mock = populate(h, ZHENG, 0.95, np.random.default_rng(1))
    exp = expected_number(h['logM'].values, ZHENG)
    assert abs(len(mock) - exp) < 5 * np.sqrt(exp)


def test_tune_density():
    h = halos()
    lm = tune_logMmin(h['logM'].values, ZHENG, 20000)
    assert np.isclose(expected_number(h['logM'].values, {**ZHENG, 'logMmin': lm}), 20000, rtol=1e-6)
    p, _ = tune_density(h['logM'].values, GHOD, 20000)
    assert np.isclose(expected_number(h['logM'].values, p), 20000, rtol=1e-6)


def test_conformity_and_fsat():
    h = halos()
    p, _ = tune_density(h['logM'].values, {**GHOD, 'conformity': True, 'fsat': 0.05}, 20000)
    mock = populate(h, p, 0.95, np.random.default_rng(2))
    assert abs(mock['is_sat'].mean() - 0.05) < 0.01


def test_nfw_inside_rvir():
    r = sample_nfw_radius(np.full(10000, 5.), np.random.default_rng(3))
    assert r.min() >= 0 and r.max() <= 1


def test_number_density_desi():
    assert 6e-4 < number_density(DATA, 'ELG1', 0.95) < 7e-4
    assert 2.5e-4 < number_density(DATA, 'LRG3', 0.95) < 3.2e-4


# ── SHAM ────────────────────────────────────────────────────────────────────
def test_standard_sham_is_rank_order():
    am = abundance_matching(objects(), 0., 500., np.random.default_rng(4))
    s = select(am, n_target=1000)
    assert len(s) == 1000 and s['Vpeak'].min() >= np.sort(am['Vpeak'].values)[-1000]
    # M_* increases with V_peak (no scatter)
    o = np.argsort(am['Vpeak'].values)
    lm = am['logMstar'].values[o]
    assert np.all(np.diff(lm[np.isfinite(lm)]) >= -1e-12)


def test_modified_sham_fsat_and_threshold():
    rng = np.random.default_rng(5)
    am = abundance_matching(objects(), 0.2, 500., rng)
    s = select(am, n_target=5000, fsat=0.1, vpeak_min=150., rng=rng)
    assert len(s) == 5000 and np.isclose(s['is_sat'].mean(), 0.1)
    assert s['Vpeakscatt'].min() > 150.


def test_weighted_draw_no_repetition():
    idx = weighted_draw(np.ones(1000), 500, np.random.default_rng(6))
    assert len(np.unique(idx)) == 500


def test_pdf_selection():
    rng = np.random.default_rng(7)
    am = abundance_matching(objects(), 0.2, 500., rng)
    s = select_pdf(am, 5000, 0.1, 150., 20., rng)
    assert np.isclose(s['is_sat'].mean(), 0.1)
    assert abs(np.median(s.loc[s['is_sat'] == 0, 'Vpeakscatt']) - 150.) < 10.   # centred on V_mean


def test_decorated_sham_tng():
    """Hydro-calibrated decoration: SMF and P(sSFR | M_*) are preserved, and with rho = 1
    the oldest halos (highest z_starve) host the quenched galaxies."""
    rng = np.random.default_rng(8)
    n = 50000
    vpeak = 80. * rng.pareto(2.5, n) + 80.
    zst = rng.uniform(0.8, 4., n)
    dark = pd.DataFrame(dict(Vpeak=vpeak, z_starve=zst, is_sat=np.zeros(n, int)))
    lms = 8.8 + 2.5 * np.log10(vpeak / 80.) + rng.normal(0, 0.1, n)
    hydro = pd.DataFrame(dict(logMstar=lms, log_sSFR=np.where(rng.uniform(size=n) < 0.4, -11.5, -9.2),
                              g_i=rng.uniform(0.3, 1.2, n)))
    dark['logMstar'] = sham_stellar_mass(dark['Vpeak'].values, hydro['logMstar'].values, 0.125, rng)
    assert np.allclose(np.sort(dark['logMstar']), np.sort(lms))           # same stellar mass function
    m = decorate(dark, hydro, 'z_starve', 1.0, (8.75, 9.3, 10.0, 12.5), rng)
    q = m['log_sSFR'] < -10
    assert abs(q.mean() - 0.4) < 0.02                                       # same quenched fraction
    assert m.loc[q, 'z_starve'].median() > m.loc[~q, 'z_starve'].median() + 1.   # assembly bias
    m0 = decorate(dark, hydro, 'z_starve', 0.0, (8.75, 9.3, 10.0, 12.5), rng)
    q0 = m0['log_sSFR'] < -10
    assert abs(m0.loc[q0, 'z_starve'].median() - m0.loc[~q0, 'z_starve'].median()) < 0.1   # rho = 0: none
    lrg = select_tracer(m, 'LRG3', 1000)
    assert len(lrg) == 1000 and (lrg['log_sSFR'] < -10).all()
