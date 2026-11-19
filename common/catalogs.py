# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
common/catalogs.py
==================
Reading the halo catalogs and writing the mocks for FCFC.

Catalogs are read through a "column map" (COLUMNS in config_common.py): the keys
are the names used by the code, the values the names in the file. Any catalog can
then be used without touching the code.

Internal names:
  x, y, z      comoving positions [h^-1 Mpc]
  vx, vy, vz   peculiar velocities [km/s]
  logM         log10 halo mass [h^-1 Msun]
  Rvir         virial radius [h^-1 Mpc, comoving]
  Vpeak        peak of the maximum circular velocity over the history [km/s]  (SHAM)
  host_id      ID of the host halo (links the particles to their host)
  conc         halo concentration                     (optional: HOD and decorated SHAM)
  is_sat       0 = host halo, 1 = satellite candidate (added by the code)
"""

from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# 1. Reading
# ─────────────────────────────────────────────────────────────────────────────
def read_halos(path, columns, boxsize=None):
    """Read an HDF5 file (one dataset per column, or a pandas table) or a
    whitespace-separated ASCII table with a header line. Returns a DataFrame
    with the internal column names.

    path     : catalog file
    columns  : dict internal_name -> name_in_the_file (None = do not read)
    boxsize  : if given, positions are wrapped into the periodic box [0, L)
    """
    path = Path(path)
    wanted = {k: v for k, v in columns.items() if v is not None}
    if path.suffix in ('.h5', '.hdf5'):
        try:
            # format of the course catalogs: one dataset per column
            with h5py.File(path, 'r') as f:
                missing = [v for v in wanted.values() if v not in f]
                if missing:
                    raise KeyError(missing)
                df = pd.DataFrame({k: f[v][:] for k, v in wanted.items()})
        except KeyError:
            # fallback: a pandas table written with DataFrame.to_hdf
            raw = pd.read_hdf(path)
            df = pd.DataFrame({k: raw[v].values for k, v in wanted.items()})
    else:
        raw = pd.read_csv(path, sep=r'\s+', comment='#')
        df = pd.DataFrame({k: raw[v].values for k, v in wanted.items()})
    if boxsize is not None:
        for c in ('x', 'y', 'z'):
            df[c] = np.mod(df[c].values, boxsize)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 2. Redshift space
# ─────────────────────────────────────────────────────────────────────────────
def zspace_factor(z, om, H0=100.):
    """a H(z) in flat LCDM. A velocity v [km/s] shifts the position by v / (a H) [h^-1 Mpc]."""
    ol = 1. - om
    return 1. / (1. + z) * H0 * np.sqrt(om * (1. + z)**3 + ol)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Writing the mock
# ─────────────────────────────────────────────────────────────────────────────
def write_mock(mock, path, z_snap, om, boxsize, fmt='hdf5'):
    """Write the mock for FCFC with redshift-space distortions along the z axis
    (plane-parallel approximation): z_s = z + v_z / (a H), wrapped into the box.
    The real-space z is kept as `zreal` (HDF5 only)."""
    x = np.mod(mock['x'].values.astype(np.float64), boxsize)
    y = np.mod(mock['y'].values.astype(np.float64), boxsize)
    zreal = mock['z'].values.astype(np.float64)
    zs = np.mod(zreal + mock['vz'].values / zspace_factor(z_snap, om), boxsize)
    path = Path(path)
    if fmt == 'hdf5':
        with h5py.File(path, 'w') as f:
            f.create_dataset('x', data=x)
            f.create_dataset('y', data=y)
            f.create_dataset('z', data=zs)                     # FCFC reads x, y, z
            f.create_dataset('zreal', data=np.mod(zreal, boxsize))
            # extra columns, useful to analyse the mock afterwards
            for c in ('vx', 'vy', 'vz', 'logM', 'is_sat', 'logMstar', 'Vpeak', 'log_sSFR'):
                if c in mock:
                    f.create_dataset(c, data=mock[c].values)
    elif fmt == 'ascii':
        np.savetxt(path, np.column_stack([x, y, zs]), fmt='%.5f', header='x y z(redshift space)')
    else:
        raise ValueError("fmt must be 'hdf5' or 'ascii'")
    return path
