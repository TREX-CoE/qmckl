"""
Python API of the QMCkl library using ctypes.

This module provides a Pythonic interface to the QMCkl (Quantum Monte Carlo
Kernel Library) C library. It uses ctypes to load the shared library at
runtime and wraps all public API functions with proper type checking and
error handling.

Usage:
    import qmckl as pq

    ctx = pq.context_create()
    pq.trexio_read(ctx, 'molecule.h5')
    mo_num = pq.get_mo_basis_mo_num(ctx)
"""

import ctypes
import ctypes.util
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Library loading
# ---------------------------------------------------------------------------

def _load_library():
    """Load the QMCkl shared library."""

    # 1. Check QMCKL_LIBDIR environment variable
    libdir = os.environ.get("QMCKL_LIBDIR", "")
    if libdir:
        for name in ("libqmckl.so", "libqmckl.dylib"):
            path = os.path.join(libdir, name)
            if os.path.isfile(path):
                return ctypes.CDLL(path)

    # 2. Try ctypes.util.find_library (searches standard paths)
    lib_path = ctypes.util.find_library("qmckl")
    if lib_path is not None:
        return ctypes.CDLL(lib_path)

    # 3. Try common names directly (relies on LD_LIBRARY_PATH / DYLD_LIBRARY_PATH)
    for name in ("libqmckl.so", "libqmckl.dylib"):
        try:
            return ctypes.CDLL(name)
        except OSError:
            continue

    raise OSError(
        "Could not find the QMCkl shared library (libqmckl.so). "
        "Make sure QMCkl is installed and either:\n"
        "  - set the QMCKL_LIBDIR environment variable, or\n"
        "  - add the library directory to LD_LIBRARY_PATH."
    )


_lib = _load_library()

# ---------------------------------------------------------------------------
# C type aliases
# ---------------------------------------------------------------------------

# qmckl_context is typedef'd as int64_t
qmckl_context = ctypes.c_int64

# qmckl_exit_code is typedef'd as int32_t
qmckl_exit_code = ctypes.c_int32
