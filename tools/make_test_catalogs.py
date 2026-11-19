# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Small synthetic catalogs with the same format as the course catalogs, to test the
code without downloading the AbacusSummit files (the clustering is NOT realistic).

    python tools/make_test_catalogs.py test_catalogs/ [box size, default 300]

Writes test_catalogs/test_halocen.h5 (host halos) and test_catalogs/test_halosat.h5
(DM particles = satellite candidates). Then, in config_common.py:
    HALO_CAT = ".../test_halocen.h5";  PARTICLE_CAT = ".../test_halosat.h5";  BOXSIZE = 300.
"""
import sys
from pathlib import Path

import h5py
import numpy as np


def write(path, cols, attrs=None):
    with h5py.File(path, 'w') as f:
        for k, v in cols.items():
            f.create_dataset(k, data=v)
        for k, v in (attrs or {}).items():
            f.attrs[k] = v
    print(f"  {path}: {len(next(iter(cols.values()))):,} rows")


def main(outdir, L=300., n_host=1_500_000, seed=1):
    rng = np.random.default_rng(seed)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # host masses from a power-law mass function, 10^11 - 10^15 h^-1 Msun
    a = -0.9
    lo, hi = 1e11, 1e15
    M = (lo**a + rng.uniform(size=n_host) * (hi**a - lo**a))**(1 / a)
    # weak clustering: hosts in clumps, more massive halos in fewer, denser clumps
    centers = rng.uniform(0, L, (20000, 3))
    which = rng.integers(0, 20000, n_host)
    pos = (centers[which] + rng.normal(0, 3., (n_host, 3)) * (M[:, None] / 1e12)**-0.1) % L
    vel = rng.normal(0, 300., (n_host, 3))
    rvir = 0.21 * (M / 1e12)**(1 / 3)                          # h^-1 Mpc, comoving
    vpeak = 160. * (M / 1e12)**0.3 * (1 + rng.normal(0, 0.05, n_host))
    conc = 8. * (M / 1e12)**-0.1 * np.exp(rng.normal(0, 0.25, n_host))
    ids = np.arange(n_host, dtype=np.int64)
    write(outdir / 'test_halocen.h5',
          dict(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], vx=vel[:, 0], vy=vel[:, 1], vz=vel[:, 2],
               logMcen=np.log10(M), Rvir=rvir, IDcen=ids, Vpeak=vpeak, Vmax=0.9 * vpeak, conc=conc),
          dict(host_fraction=1.0))

    # particles: a few per halo, more in massive halos, inside R_vir, with a sampled V_peak
    npart = rng.poisson(3 + 20 * (M / 1e13)**0.9)
    h = np.repeat(np.arange(n_host), npart)
    r = rng.uniform(0, 1, len(h))**0.5 * rvir[h]
    d = rng.normal(size=(len(h), 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    ppos = (pos[h] + d * r[:, None]) % L
    pvel = vel[h] + rng.normal(0, 1., (len(h), 3)) * (300. * (M[h] / 1e13)**(1 / 3))[:, None]
    pvpeak = vpeak[h] * rng.uniform(0.2, 0.7, len(h))
    write(outdir / 'test_halosat.h5',
          dict(x=ppos[:, 0], y=ppos[:, 1], z=ppos[:, 2], vx=pvel[:, 0], vy=pvel[:, 1], vz=pvel[:, 2],
               IDcen=ids[h], Vpeak=pvpeak))


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 300.)
