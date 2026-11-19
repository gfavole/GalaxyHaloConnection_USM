# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
hod/model.py
============
HOD (Halo Occupation Distribution) models: how many central and satellite galaxies
a halo of mass M hosts on average, and where they are placed.

Three shapes of the central occupation (p['model']):

'zheng'  standard 5-parameter HOD (Zheng et al. 2007, arXiv:astro-ph/0703457):

    <N_cen>(M) = 1/2 [1 + erf((log M - log M_min) / sigma_logM)]
    <N_sat>(M) = <N_cen>(M) ((M - M_0) / M_1)^alpha        for M > M_0, else 0

    It assumes a COMPLETE sample: every massive enough halo hosts a central.

'ghod'   symmetric Gaussian for ELGs (Rocher et al. 2023, arXiv:2306.06319, Eq. 3.1):

    <N_cen>(M) = A_c / (sqrt(2 pi) sigma_M) exp[-(log M - log M_c)^2 / (2 sigma_M^2)]
    <N_sat>(M) = A_s ((M - M_0) / M_1)^alpha                for M > M_0, else 0

    INCOMPLETE sample: ELG centrals live in a range of halo masses and disappear in
    massive halos, where centrals are quenched.

'mhmq'   asymmetric Gaussian (modified High Mass Quenched, Rocher+23, Eq. 3.4):

    <N_cen>(M) = GHOD(M) [1 + erf(gamma (log M - log M_c) / (sqrt(2) sigma_M))]

    with a tail towards high masses (gamma > 0). gamma = 0 gives back the GHOD (times 2).

Satellite options (keys of p, any model):
  conformity   True = strict conformity: a halo gets satellites only if it has a central
  f_sigma_v    velocity bias: satellite velocity dispersion = f_sigma_v sigma_v(halo)
  f_exp, tau, lambda_nfw   modified profile of Rocher+23: a fraction f_exp of the
               satellites follows dN/dr ~ exp(-r / (tau r_s)) (they can sit outside
               R_vir); the rest an NFW with r_s -> r_s / lambda_nfw (less concentrated
               if lambda_nfw < 1)
  placement    'nfw' forces the analytic profile for this model (needed for the
               velocity bias and the modified profile), whatever SAT_PLACEMENT says

Two ways of placing the satellites (SAT_PLACEMENT):
  'nfw'        analytic NFW profile + Gaussian velocities sigma_1D = V_vir / sqrt(2)
  'particles'  on dark-matter particles of the same halo (linked by host_id), with
               their positions and velocities. If a halo has fewer particles than
               satellites, all its particles are used (and the shortfall is reported).
"""

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import erf

G_NEWTON = 4.30091e-9          # Mpc (km/s)^2 / Msun; with M in h^-1 Msun and R in h^-1 Mpc the h cancel


# ─────────────────────────────────────────────────────────────────────────────
# 1. Mean occupations
# ─────────────────────────────────────────────────────────────────────────────
def ncen_mean(logM, logMmin, sigma_logM):
    """Zheng centrals: a smoothed step function that saturates at 1."""
    return 0.5 * (1. + erf((logM - logMmin) / sigma_logM))


def nsat_mean(logM, logMmin, sigma_logM, logM0, logM1, alpha, times_ncen=True):
    """Zheng satellites: power law, multiplied by <N_cen> (times_ncen=True)."""
    M = 10**np.asarray(logM, dtype=float)
    x = np.clip((M - 10**logM0) / 10**logM1, 0., None)   # 0 for M <= M_0
    ns = x**alpha
    return ns * ncen_mean(logM, logMmin, sigma_logM) if times_ncen else ns


def ncen_ghod(logM, Ac, logMc, sigma_M):
    """Gaussian centrals in log M (Rocher+23): peak at M_c, width sigma_M, amplitude A_c."""
    return Ac / (np.sqrt(2. * np.pi) * sigma_M) * np.exp(-(logM - logMc)**2 / (2. * sigma_M**2))


def nsat_powerlaw(logM, As, logM0, logM1, alpha):
    """Rocher+23 satellites: power law, independent of the centrals."""
    x = (10**np.asarray(logM, dtype=float) - 10**logM0) / 10**logM1
    out = np.zeros_like(x)
    pos = x > 0                                           # only halos with M > M_0
    out[pos] = As * x[pos]**alpha                         # (alpha can be negative)
    return out


def mean_occupations(logM, p):
    """(<N_cen>, <N_sat>) per halo. With conformity, <N_sat> is the mean number of
    satellites in a halo that DOES have a central."""
    model = p.get('model', 'zheng')
    if model == 'zheng':
        nc = ncen_mean(logM, p['logMmin'], p['sigma_logM'])
        ns = nsat_mean(logM, p['logMmin'], p['sigma_logM'], p['logM0'], p['logM1'], p['alpha'],
                       p.get('times_ncen', True))
    elif model in ('ghod', 'mhmq'):
        nc = ncen_ghod(logM, p['Ac'], p['logMc'], p['sigma_M'])
        if model == 'mhmq':                       # asymmetry term
            nc = nc * (1. + erf(p['gamma'] * (logM - p['logMc']) / (np.sqrt(2.) * p['sigma_M'])))
        ns = nsat_powerlaw(logM, p['As'], p['logM0'], p['logM1'], p['alpha'])
    else:
        raise ValueError("p['model'] must be 'zheng', 'ghod' or 'mhmq'")
    return np.clip(nc, 0., 1.), ns                # a probability cannot exceed 1


def expected_number(logM, p):
    """Expected total number of galaxies for a halo catalog (array of log M)."""
    nc, ns = mean_occupations(logM, p)
    return nc.sum() + ((nc * ns).sum() if p.get('conformity', False) else ns.sum())


# ─────────────────────────────────────────────────────────────────────────────
# 2. Matching the DESI number density
# ─────────────────────────────────────────────────────────────────────────────
def tune_logMmin(logM, p, n_target):
    """Zheng: log M_min such that the expected number of galaxies is n_target."""
    f = lambda lm: expected_number(logM, {**p, 'logMmin': lm}) - n_target
    return brentq(f, 10., 16.)


def tune_amplitudes(logM, p, n_target):
    """GHOD/mHMQ: a common factor k on A_c and A_s (as in Rocher+23). It changes the
    density without changing the satellite fraction (without conformity)."""
    f = lambda k: expected_number(logM, {**p, 'Ac': k * p['Ac'], 'As': k * p['As']}) - n_target
    k = brentq(f, 1e-8, 1e4)
    nc, _ = mean_occupations(logM, {**p, 'Ac': k * p['Ac']})
    if np.any(nc >= 1.):
        print("    warning: <N_cen> reaches 1 in some halos (clipped): the density is hard to reach")
    return {**p, 'Ac': k * p['Ac'], 'As': k * p['As']}, k


def tune_amplitudes_fsat(logM, p, n_target, fsat, n_iter=20):
    """GHOD/mHMQ: A_c and A_s set separately to give n_target galaxies AND a satellite
    fraction fsat (iterative, because with conformity the satellites depend on A_c)."""
    p = dict(p)
    for _ in range(n_iter):
        nc, ns = mean_occupations(logM, p)
        ncen = nc.sum()
        nsat = (nc * ns).sum() if p.get('conformity', False) else ns.sum()
        p['Ac'] *= (1. - fsat) * n_target / ncen
        p['As'] *= fsat * n_target / max(nsat, 1e-30)
    return p


def tune_density(logM, p, n_target):
    """Match the DESI density: log M_min (Zheng) or the amplitudes A_c, A_s (GHOD/mHMQ).
    If p['fsat'] is given, the satellite fraction is fixed too. Returns (p, message)."""
    if p.get('model', 'zheng') == 'zheng':
        lm = tune_logMmin(logM, p, n_target)
        return {**p, 'logMmin': lm}, f"log M_min = {lm:.3f}"
    if p.get('fsat') is not None:
        p2 = tune_amplitudes_fsat(logM, p, n_target, p['fsat'])
        return p2, f"A_c = {p2['Ac']:.4f}, A_s = {p2['As']:.4f} (target f_sat = {p['fsat']:.3f})"
    p2, k = tune_amplitudes(logM, p, n_target)
    return p2, f"A_c = {p2['Ac']:.4f}, A_s = {p2['As']:.4f} (factor {k:.4f})"


# ─────────────────────────────────────────────────────────────────────────────
# 3. Satellite profile
# ─────────────────────────────────────────────────────────────────────────────
def concentration_duffy08(logM, z):
    """c_vir(M, z) from Duffy et al. 2008 (full sample, virial definition)."""
    return 7.85 * (10**logM / 2e12)**(-0.081) * (1. + z)**(-0.71)


def sample_nfw_radius(c, rng):
    """r / R_vir drawn from an NFW profile truncated at R_vir (inverse CDF on a grid).
    The enclosed mass is m(y) = ln(1 + y) - y / (1 + y), with y = c r / R_vir."""
    c = np.atleast_1d(c)
    out = np.empty(len(c))
    m = lambda y: np.log1p(y) - y / (1. + y)
    u = rng.uniform(size=len(c))
    grid = np.linspace(0., 1., 513)[1:]
    for cc in np.unique(np.round(c, 2)):          # one CDF per concentration value
        sel = np.round(c, 2) == cc
        cdf = m(cc * grid) / m(cc)
        out[sel] = np.interp(u[sel], np.r_[0., cdf], np.r_[0., grid])
    return out


def random_directions(n, rng):
    """n isotropic unit vectors."""
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1)[:, None]


def _satellites_nfw(h, z_snap, rng, f_sigma_v=1., f_exp=0., tau=None, lambda_nfw=1.):
    """Satellites on an NFW profile with Gaussian velocities. h has one row per
    satellite (that of its host). Optional: velocity bias (f_sigma_v) and the
    modified profile of Rocher+23 (f_exp, tau, lambda_nfw)."""
    R = h['Rvir'].values
    c = h['conc'].values if 'conc' in h else concentration_duffy08(h['logM'].values, z_snap)
    # radius: NFW with r_s -> r_s / lambda_nfw
    r = sample_nfw_radius(c, rng) * R / lambda_nfw
    if f_exp > 0:                                  # a fraction f_exp on an exponential profile
        expo = rng.uniform(size=len(h)) < f_exp
        r[expo] = rng.exponential(tau * R[expo] / c[expo])
    d = random_directions(len(h), rng) * r[:, None]
    # velocity: that of the halo + Gaussian dispersion sigma_1D = f_sigma_v V_vir / sqrt(2)
    a = 1. / (1. + z_snap)
    vvir = np.sqrt(G_NEWTON * 10**h['logM'].values / (R * a))       # physical R = a x comoving R
    dv = rng.normal(size=(len(h), 3)) * (f_sigma_v * vvir / np.sqrt(2.))[:, None]
    return pd.DataFrame({'x': h['x'].values + d[:, 0], 'y': h['y'].values + d[:, 1],
                         'z': h['z'].values + d[:, 2],
                         'vx': h['vx'].values + dv[:, 0], 'vy': h['vy'].values + dv[:, 1],
                         'vz': h['vz'].values + dv[:, 2], 'logM': h['logM'].values, 'is_sat': 1})


def _satellites_particles(host_ids, nsat, particles, logM, rng):
    """nsat[i] particles (without repetition) among those of host host_ids[i].

    Returns the placed satellites and, per host, how many could not be placed
    because the host has fewer particles than satellites (an array aligned
    with host_ids). populate() puts those on an NFW profile."""
    sel = nsat > 0
    want = pd.Series(nsat[sel], index=host_ids[sel])
    want_logM = pd.Series(logM[sel], index=host_ids[sel])
    part = particles[particles['host_id'].isin(want.index)]
    # random order within each host, then keep the first nsat
    part = part.iloc[rng.permutation(len(part))]
    rank = part.groupby('host_id').cumcount().values
    keep = rank < want.reindex(part['host_id']).values
    out = part.loc[keep, ['x', 'y', 'z', 'vx', 'vy', 'vz']].copy()
    out['logM'] = want_logM.reindex(part.loc[keep, 'host_id']).values
    out['is_sat'] = 1
    # satellites still missing in each host
    placed = part.loc[keep, 'host_id'].value_counts()
    missing = (want - placed.reindex(want.index).fillna(0).astype(int)).clip(lower=0)
    n_missing = np.zeros(len(host_ids), dtype=int)
    n_missing[sel] = missing.values
    return out.reset_index(drop=True), n_missing


# ─────────────────────────────────────────────────────────────────────────────
# 4. Populating the halos
# ─────────────────────────────────────────────────────────────────────────────
def populate(halos, p, z_snap, rng, placement='nfw', particles=None):
    """Build a galaxy mock from a host-halo catalog.

    halos     : DataFrame with x, y, z, vx, vy, vz, logM, Rvir [h^-1 Mpc] (optional conc, host_id)
    p         : model parameters (see above)
    placement : 'nfw' or 'particles' (then `particles` with host_id, x, y, z, vx, vy, vz)
    """
    if 'is_sat' in halos:
        halos = halos[halos['is_sat'] == 0]
    logM = halos['logM'].values

    # centrals: one Bernoulli draw per halo with probability <N_cen>
    pc, lam = mean_occupations(logM, p)
    is_cen = rng.uniform(size=len(halos)) < pc
    cen = halos.loc[is_cen, ['x', 'y', 'z', 'vx', 'vy', 'vz', 'logM']].copy()
    cen['is_sat'] = 0

    # satellites: Poisson number with mean <N_sat>
    nsat = rng.poisson(lam)
    if p.get('conformity', False):         # strict conformity: no central, no satellites
        nsat[~is_cen] = 0
    placement = p.get('placement', placement)
    if placement == 'nfw':
        sat = _satellites_nfw(halos.iloc[np.repeat(np.arange(len(halos)), nsat)], z_snap, rng,
                              f_sigma_v=p.get('f_sigma_v', 1.), f_exp=p.get('f_exp', 0.),
                              tau=p.get('tau'), lambda_nfw=p.get('lambda_nfw', 1.))
    elif placement == 'particles':
        sat, n_missing = _satellites_particles(halos['host_id'].values, nsat, particles, logM, rng)
        if n_missing.sum() > 0:
            # hosts with fewer particles than satellites: the rest go on an NFW profile
            print(f"    note: {int(n_missing.sum()):,} of {int(nsat.sum()):,} satellites placed on NFW "
                  f"(hosts with fewer particles than satellites)")
            extra = _satellites_nfw(halos.iloc[np.repeat(np.arange(len(halos)), n_missing)], z_snap, rng,
                                    f_sigma_v=p.get('f_sigma_v', 1.), f_exp=p.get('f_exp', 0.),
                                    tau=p.get('tau'), lambda_nfw=p.get('lambda_nfw', 1.))
            sat = pd.concat([sat, extra], ignore_index=True)
    else:
        raise ValueError("placement must be 'nfw' or 'particles'")
    return pd.concat([cen, sat], ignore_index=True)
