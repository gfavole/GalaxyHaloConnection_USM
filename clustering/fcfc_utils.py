# Copyright (c) 2026 Ginevra Favole <gifavole@ull.edu.es>
# Instituto de Astrofisica de Canarias (IAC) and Universidad de La Laguna (ULL)
# Course "Connecting galaxies to their dark matter halos", USM, November 2026
#
"""
clustering/fcfc_utils.py
========================
Interface to FCFC (C. Zhao, https://github.com/cheng-zhao/FCFC), the C code that
counts pairs and computes the correlation functions in a periodic box.

How it works: the .conf templates in clustering/fcfc_inputs/ are copied to a
temporary file in which only the keys that change at each run are rewritten
(input catalog, output files, box size, binning file, catalog format); then
./FCFC_2PT_BOX is run from FCFC_DIR.

Binning (the same as the DESI DR1 measurements in data/desi_dr1):
  26 log bins in s (and rp) in [0.17, 32] h^-1 Mpc   (fcfc_inputs/logbins.txt)
  120 mu bins; pi in [0, 40] h^-1 Mpc with d pi = 2 for wp
  natural estimator DD/RR - 1 with analytic RR in a periodic box.

The mock must already be in redshift space along z (see common/catalogs.py).
"""

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import interpolate

FCFC_INPUTS = Path(__file__).resolve().parent / 'fcfc_inputs'
SEP_BIN_FILE = FCFC_INPUTS / 'logbins.txt'
# one template per correlation type (auto or cross) and statistic (multipoles or wp)
CONF = {('auto', 'xil'): FCFC_INPUTS / 'fcfc_auto.conf', ('auto', 'wp'): FCFC_INPUTS / 'fcfc_auto_wp.conf',
        ('cross', 'xil'): FCFC_INPUTS / 'fcfc_cross.conf', ('cross', 'wp'): FCFC_INPUTS / 'fcfc_cross_wp.conf'}


def _format_keys(fmt, n):
    """FCFC keys that depend on the catalog format, for n catalogs.
    'hdf5': columns by name (FCFC compiled with WITH_HDF5 := T); 'ascii': by position."""
    if fmt == 'hdf5':
        return {'CATALOG_TYPE': '[' + ','.join(['2'] * n) + ']',
                'POSITION': '[' + ', '.join(['${x}, ${y}, ${z}'] * n) + ']'}
    if fmt == 'ascii':
        return {'CATALOG_TYPE': '[' + ','.join(['0'] * n) + ']',
                'POSITION': '[' + ', '.join(['$1, $2, $3'] * n) + ']',
                'ASCII_SKIP': '[' + ','.join(['0'] * n) + ']',
                'ASCII_COMMENT': '[' + ','.join(["'#'"] * n) + ']',
                'ASCII_FORMATTER': '[' + ','.join(['"%lf %lf %lf"'] * n) + ']'}
    raise ValueError("fmt must be 'hdf5' or 'ascii'")


def _format_fcfc_value(val):
    """Python lists -> FCFC syntax '[ a, b ]'."""
    if isinstance(val, (list, tuple)):
        return "[ " + ", ".join(str(v) for v in val) + " ]"
    return str(val)


def _run(base_conf_file, update_map, fcfc_dir, nthreads=4):
    """Write a temporary .conf from the template and run FCFC."""
    lines = []
    with open(base_conf_file) as f:
        for line in f:
            stripped = line.strip()
            if '=' in stripped and not stripped.startswith('#'):
                key = stripped.split('=')[0].strip()
                if key in update_map:                     # key that changes in this run
                    lines.append(f"{line.split('=')[0]}= {_format_fcfc_value(update_map[key])}\n")
                    continue
            lines.append(line)
    fd, tmp = tempfile.mkstemp(suffix='.conf')
    with os.fdopen(fd, 'w') as f:
        f.writelines(lines)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(nthreads)                # OpenMP threads
    exe = Path(fcfc_dir) / 'FCFC_2PT_BOX'
    if not exe.exists():
        raise FileNotFoundError(f"{exe} not found: check FCFC_DIR in config_common.py")
    try:
        subprocess.run(["./FCFC_2PT_BOX", "--conf", tmp], cwd=fcfc_dir, check=True,
                       capture_output=True, text=True, env=env)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"FCFC failed:\n{e.stdout}\n{e.stderr}") from e
    finally:
        os.remove(tmp)


def measure(mock_files, outdir, fcfc_dir, box_size, scale=None, fmt='hdf5', nthreads=4, tag='mock'):
    """Measure xi0, xi2 (s) and wp (rp) with FCFC in a periodic box.

    mock_files : [mock] for an auto-correlation, [ELG, LRG] for the cross (estimator DS/RR - 1)
    scale      : if given, everything is interpolated at these separations (those of the data)
    Returns (s, xi0, xi2, rp, wp)."""
    mock_files = [str(m) for m in np.atleast_1d(mock_files)]
    kind = 'auto' if len(mock_files) == 1 else 'cross'
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pc = ([str(outdir / f'{tag}_{x}.txt') for x in ('DD', 'DS', 'SS')] if kind == 'cross'
          else [str(outdir / f'{tag}_DD.txt')])
    cf = str(outdir / f'{tag}_cf.txt')
    xil, wpf = str(outdir / f'{tag}_xil.txt'), str(outdir / f'{tag}_wp.txt')
    keys = {'CATALOG': mock_files, 'BOX_SIZE': float(box_size), 'SEP_BIN_FILE': str(SEP_BIN_FILE),
            'PAIR_COUNT_FILE': pc, 'CF_OUTPUT_FILE': cf, **_format_keys(fmt, len(mock_files))}

    def _clean():                                  # remove the intermediate pair-count files
        for f in pc + [cf]:
            if os.path.exists(f):
                os.remove(f)

    # multipoles xi0, xi2
    _run(CONF[(kind, 'xil')], {**keys, 'MULTIPOLE_FILE': xil}, fcfc_dir, nthreads)
    s, _, _, xi0, xi2, _ = np.loadtxt(xil, comments='#').T
    _clean()
    # projected correlation function wp
    _run(CONF[(kind, 'wp')], {**keys, 'PROJECTED_FILE': wpf}, fcfc_dir, nthreads)
    rp, _, _, wp = np.loadtxt(wpf, comments='#').T
    _clean()

    if scale is None:
        return s, xi0, xi2, rp, wp
    itp = lambda x, y: interpolate.interp1d(x, y, fill_value="extrapolate")(scale)
    return scale, itp(s, xi0), itp(s, xi2), scale, itp(rp, wp)
