# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
Fake TNG inputs with the same format as the TNG300 catalogs of the course Drive, to
test 05_SHAM_decorated without them (NOT realistic):

    python tools/tng/make_fake_tng_inputs.py test_catalogs/
"""
import sys
from pathlib import Path

import h5py
import numpy as np


def write(path, cols, attrs):
    with h5py.File(path, 'w') as f:
        for k, v in cols.items():
            f.create_dataset(k, data=v)
        for k, v in attrs.items():
            f.attrs[k] = v
    print(f"  {path}: {len(next(iter(cols.values()))):,} rows")


def main(outdir, L=205., n=300_000, seed=2):
    rng = np.random.default_rng(seed)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    attrs = dict(boxsize=L, redshift=0.82, h=0.6774, Om=0.3089, snap=55)
    # dark subhalos in clumps; older (higher z_starve) subhalos in denser clumps
    centers = rng.uniform(0, L, (8000, 3))
    which = rng.integers(0, 8000, n)
    pos = (centers[which] + rng.normal(0, 2., (n, 3))) % L
    vpeak = 80. * rng.pareto(2.5, n) + 80.
    zst = np.clip(0.82 + 1.5 * (which / 8000.) + rng.normal(0, 0.4, n), 0.82, None)
    is_sat = (rng.uniform(size=n) < 0.2).astype(np.int8)
    vel = rng.normal(0, 300., (n, 3))
    write(out / 'TNG300-1-Dark_snap055_subhalos.h5',
          dict(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], vx=vel[:, 0], vy=vel[:, 1], vz=vel[:, 2],
               is_sat=is_sat, host_id=which.astype(np.int64), logMsub=10 + 3 * np.log10(vpeak / 80.),
               Vmax=vpeak * 0.9, Vpeak=vpeak, z_half=zst, z_char=np.full(n, 0.82), z_acc=np.full(n, 0.82),
               z_starve=zst, c_infall=1.2 + 0.1 * zst + rng.normal(0, 0.1, n)), attrs)
    # hydro galaxies: the same objects, M_* from V_peak, quenched if old
    lms = 8.6 + 2.6 * np.log10(vpeak / 80.) + rng.normal(0, 0.15, n)
    q = zst + rng.normal(0, 0.3, n) > 1.6
    ssfr = np.where(q, -11.5 + 0.4 * rng.normal(size=n), -9.2 + 0.3 * rng.normal(size=n))
    keep = lms > 8.75
    write(out / 'TNG300-1_snap055_galaxies.h5',
          dict(x=pos[keep, 0], y=pos[keep, 1], z=pos[keep, 2], vx=vel[keep, 0], vy=vel[keep, 1], vz=vel[keep, 2],
               is_sat=is_sat[keep], host_id=which[keep].astype(np.int64), logMstar=lms[keep],
               SFR=10**(ssfr[keep] + lms[keep]), log_sSFR=ssfr[keep], g_i=np.where(q, 1.1, 0.5)[keep]), attrs)


if __name__ == '__main__':
    main(sys.argv[1])
