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

# ---------------------------------------------------------------------------
# Error codes
# ---------------------------------------------------------------------------

QMCKL_SUCCESS              = 0
QMCKL_INVALID_ARG_1        = 1
QMCKL_INVALID_ARG_2        = 2
QMCKL_INVALID_ARG_3        = 3
QMCKL_INVALID_ARG_4        = 4
QMCKL_INVALID_ARG_5        = 5
QMCKL_INVALID_ARG_6        = 6
QMCKL_INVALID_ARG_7        = 7
QMCKL_INVALID_ARG_8        = 8
QMCKL_INVALID_ARG_9        = 9
QMCKL_INVALID_ARG_10       = 10
QMCKL_INVALID_ARG_11       = 11
QMCKL_INVALID_ARG_12       = 12
QMCKL_INVALID_ARG_13       = 13
QMCKL_INVALID_ARG_14       = 14
QMCKL_INVALID_ARG_15       = 15
QMCKL_INVALID_ARG_16       = 16
QMCKL_INVALID_ARG_17       = 17
QMCKL_INVALID_ARG_18       = 18
QMCKL_INVALID_ARG_19       = 19
QMCKL_INVALID_ARG_20       = 20
QMCKL_FAILURE              = 101
QMCKL_ERRNO                = 102
QMCKL_INVALID_CONTEXT      = 103
QMCKL_ALLOCATION_FAILED    = 104
QMCKL_DEALLOCATION_FAILED  = 105
QMCKL_NOT_PROVIDED         = 106
QMCKL_OUT_OF_BOUNDS        = 107
QMCKL_ALREADY_SET          = 108
QMCKL_INVALID_EXIT_CODE    = 109

# ---------------------------------------------------------------------------
# Low-level C function signatures
# ---------------------------------------------------------------------------

# const char* qmckl_string_of_error(const qmckl_exit_code error);
_lib.qmckl_string_of_error.argtypes = [qmckl_exit_code]
_lib.qmckl_string_of_error.restype  = ctypes.c_char_p

# --- Context ---

# qmckl_context qmckl_context_create();
_lib.qmckl_context_create.argtypes = []
_lib.qmckl_context_create.restype  = qmckl_context

# qmckl_exit_code qmckl_context_destroy(const qmckl_context context);
_lib.qmckl_context_destroy.argtypes = [qmckl_context]
_lib.qmckl_context_destroy.restype  = qmckl_exit_code

# qmckl_context qmckl_context_copy(const qmckl_context context);
_lib.qmckl_context_copy.argtypes = [qmckl_context]
_lib.qmckl_context_copy.restype  = qmckl_context

# --- TREXIO ---

# qmckl_exit_code qmckl_trexio_read(const qmckl_context context,
#                                    const char* file_name,
#                                    const int64_t size_max);
_lib.qmckl_trexio_read.argtypes = [qmckl_context, ctypes.c_char_p, ctypes.c_int64]
_lib.qmckl_trexio_read.restype  = qmckl_exit_code

# --- Numerical precision ---

# qmckl_exit_code qmckl_set_numprec_precision(const qmckl_context context, const int precision);
_lib.qmckl_set_numprec_precision.argtypes = [qmckl_context, ctypes.c_int32]
_lib.qmckl_set_numprec_precision.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_numprec_range(const qmckl_context context, const int range);
_lib.qmckl_set_numprec_range.argtypes = [qmckl_context, ctypes.c_int32]
_lib.qmckl_set_numprec_range.restype  = qmckl_exit_code

# --- Electron ---

# qmckl_exit_code qmckl_set_electron_num(qmckl_context, int64_t, int64_t);
_lib.qmckl_set_electron_num.argtypes = [qmckl_context, ctypes.c_int64, ctypes.c_int64]
_lib.qmckl_set_electron_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_electron_coord(qmckl_context, char, int64_t, const double*, int64_t);
_lib.qmckl_set_electron_coord.argtypes = [
    qmckl_context, ctypes.c_char, ctypes.c_int64,
    ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_electron_coord.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_num(const qmckl_context, int64_t*);
_lib.qmckl_get_electron_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_electron_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_up_num(const qmckl_context, int64_t*);
_lib.qmckl_get_electron_up_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_electron_up_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_down_num(const qmckl_context, int64_t*);
_lib.qmckl_get_electron_down_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_electron_down_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_walk_num(const qmckl_context, int64_t*);
_lib.qmckl_get_electron_walk_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_electron_walk_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_coord(const qmckl_context, char, double*, int64_t);
_lib.qmckl_get_electron_coord.argtypes = [
    qmckl_context, ctypes.c_char,
    ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_electron_coord.restype  = qmckl_exit_code

# --- Nucleus ---

# qmckl_exit_code qmckl_get_nucleus_num(const qmckl_context, int64_t*);
_lib.qmckl_get_nucleus_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_nucleus_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_nucleus_charge(const qmckl_context, double*, int64_t);
_lib.qmckl_get_nucleus_charge.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_nucleus_charge.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_nucleus_coord(const qmckl_context, char, double*, int64_t);
_lib.qmckl_get_nucleus_coord.argtypes = [
    qmckl_context, ctypes.c_char,
    ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_nucleus_coord.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_nucleus_num(qmckl_context, int64_t);
_lib.qmckl_set_nucleus_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_nucleus_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_nucleus_charge(qmckl_context, const double*, int64_t);
_lib.qmckl_set_nucleus_charge.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_nucleus_charge.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_nucleus_coord(qmckl_context, char, const double*, int64_t);
_lib.qmckl_set_nucleus_coord.argtypes = [
    qmckl_context, ctypes.c_char,
    ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_nucleus_coord.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_nucleus_nn_distance(qmckl_context, double*);
_lib.qmckl_get_nucleus_nn_distance.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double)
]
_lib.qmckl_get_nucleus_nn_distance.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_nucleus_repulsion(qmckl_context, double*);
_lib.qmckl_get_nucleus_repulsion.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double)
]
_lib.qmckl_get_nucleus_repulsion.restype  = qmckl_exit_code

# --- Point ---

# qmckl_exit_code qmckl_get_point_num(const qmckl_context, int64_t*);
_lib.qmckl_get_point_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_point_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_point(const qmckl_context, char, double*, int64_t);
_lib.qmckl_get_point.argtypes = [
    qmckl_context, ctypes.c_char,
    ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_point.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_point(qmckl_context, char, int64_t, const double*, int64_t);
_lib.qmckl_set_point.argtypes = [
    qmckl_context, ctypes.c_char, ctypes.c_int64,
    ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_point.restype  = qmckl_exit_code

# --- AO basis ---

# qmckl_exit_code qmckl_get_ao_basis_type(const qmckl_context, char*);
_lib.qmckl_get_ao_basis_type.argtypes = [qmckl_context, ctypes.c_char_p]
_lib.qmckl_get_ao_basis_type.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_type(qmckl_context, const char*, int64_t);
_lib.qmckl_set_ao_basis_type.argtypes = [qmckl_context, ctypes.c_char_p, ctypes.c_int64]
_lib.qmckl_set_ao_basis_type.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_shell_num(qmckl_context, int64_t);
_lib.qmckl_set_ao_basis_shell_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_ao_basis_shell_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_prim_num(qmckl_context, int64_t);
_lib.qmckl_set_ao_basis_prim_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_ao_basis_prim_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_ao_num(qmckl_context, int64_t);
_lib.qmckl_set_ao_basis_ao_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_ao_basis_ao_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_nucleus_index(qmckl_context, const int64_t*, int64_t);
_lib.qmckl_set_ao_basis_nucleus_index.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_int64), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_nucleus_index.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_nucleus_shell_num(qmckl_context, const int64_t*, int64_t);
_lib.qmckl_set_ao_basis_nucleus_shell_num.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_int64), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_nucleus_shell_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_shell_ang_mom(qmckl_context, const int32_t*, int64_t);
_lib.qmckl_set_ao_basis_shell_ang_mom.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_int32), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_shell_ang_mom.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_shell_prim_num(qmckl_context, const int64_t*, int64_t);
_lib.qmckl_set_ao_basis_shell_prim_num.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_int64), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_shell_prim_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_shell_prim_index(qmckl_context, const int64_t*, int64_t);
_lib.qmckl_set_ao_basis_shell_prim_index.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_int64), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_shell_prim_index.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_shell_factor(qmckl_context, const double*, int64_t);
_lib.qmckl_set_ao_basis_shell_factor.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_shell_factor.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_exponent(qmckl_context, const double*, int64_t);
_lib.qmckl_set_ao_basis_exponent.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_exponent.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_coefficient(qmckl_context, const double*, int64_t);
_lib.qmckl_set_ao_basis_coefficient.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_coefficient.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_prim_factor(qmckl_context, const double*, int64_t);
_lib.qmckl_set_ao_basis_prim_factor.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_prim_factor.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_ao_basis_ao_factor(qmckl_context, const double*, int64_t);
_lib.qmckl_set_ao_basis_ao_factor.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_ao_basis_ao_factor.restype  = qmckl_exit_code

# --- MO basis ---

# qmckl_exit_code qmckl_get_mo_basis_mo_num(const qmckl_context, int64_t*);
_lib.qmckl_get_mo_basis_mo_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_mo_basis_mo_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_mo_basis_mo_num(qmckl_context, int64_t);
_lib.qmckl_set_mo_basis_mo_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_mo_basis_mo_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_mo_basis_coefficient(const qmckl_context, double*, int64_t);
_lib.qmckl_get_mo_basis_coefficient.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_mo_basis_coefficient.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_mo_basis_coefficient(qmckl_context, const double*, int64_t);
_lib.qmckl_set_mo_basis_coefficient.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_mo_basis_coefficient.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_mo_basis_mo_vgl(qmckl_context, double*, int64_t);
_lib.qmckl_get_mo_basis_mo_vgl.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_mo_basis_mo_vgl.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_mo_basis_mo_vgl_inplace(qmckl_context, double*, int64_t);
_lib.qmckl_get_mo_basis_mo_vgl_inplace.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_get_mo_basis_mo_vgl_inplace.restype  = qmckl_exit_code

# --- Jastrow CHAMP ---

# qmckl_exit_code qmckl_set_jastrow_champ_rescale_factor_ee(qmckl_context, double);
_lib.qmckl_set_jastrow_champ_rescale_factor_ee.argtypes = [qmckl_context, ctypes.c_double]
_lib.qmckl_set_jastrow_champ_rescale_factor_ee.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_jastrow_champ_rescale_factor_en(qmckl_context, const double*, int64_t);
_lib.qmckl_set_jastrow_champ_rescale_factor_en.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double), ctypes.c_int64
]
_lib.qmckl_set_jastrow_champ_rescale_factor_en.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_jastrow_champ_aord_num(qmckl_context, int64_t);
_lib.qmckl_set_jastrow_champ_aord_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_jastrow_champ_aord_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_jastrow_champ_bord_num(qmckl_context, int64_t);
_lib.qmckl_set_jastrow_champ_bord_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_jastrow_champ_bord_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_set_jastrow_champ_cord_num(qmckl_context, int64_t);
_lib.qmckl_set_jastrow_champ_cord_num.argtypes = [qmckl_context, ctypes.c_int64]
_lib.qmckl_set_jastrow_champ_cord_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_jastrow_champ_aord_num(qmckl_context, int64_t*);
_lib.qmckl_get_jastrow_champ_aord_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_jastrow_champ_aord_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_jastrow_champ_bord_num(qmckl_context, int64_t*);
_lib.qmckl_get_jastrow_champ_bord_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_jastrow_champ_bord_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_jastrow_champ_cord_num(qmckl_context, int64_t*);
_lib.qmckl_get_jastrow_champ_cord_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_jastrow_champ_cord_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_jastrow_champ_type_nucl_num(qmckl_context, int64_t*);
_lib.qmckl_get_jastrow_champ_type_nucl_num.argtypes = [qmckl_context, ctypes.POINTER(ctypes.c_int64)]
_lib.qmckl_get_jastrow_champ_type_nucl_num.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_jastrow_champ_type_nucl_vector(qmckl_context, int64_t*, int64_t);
_lib.qmckl_get_jastrow_champ_type_nucl_vector.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_int64), ctypes.c_int64
]
_lib.qmckl_get_jastrow_champ_type_nucl_vector.restype  = qmckl_exit_code

# --- Distance ---

# qmckl_exit_code qmckl_get_electron_ee_distance(qmckl_context, double*);
_lib.qmckl_get_electron_ee_distance.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double)
]
_lib.qmckl_get_electron_ee_distance.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_en_distance(qmckl_context, double*);
_lib.qmckl_get_electron_en_distance.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double)
]
_lib.qmckl_get_electron_en_distance.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_ee_potential(qmckl_context, double*);
_lib.qmckl_get_electron_ee_potential.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double)
]
_lib.qmckl_get_electron_ee_potential.restype  = qmckl_exit_code

# qmckl_exit_code qmckl_get_electron_en_potential(qmckl_context, double*);
_lib.qmckl_get_electron_en_potential.argtypes = [
    qmckl_context, ctypes.POINTER(ctypes.c_double)
]
_lib.qmckl_get_electron_en_potential.restype  = qmckl_exit_code


# ---------------------------------------------------------------------------
# Error handling helper
# ---------------------------------------------------------------------------

def _check_rc(rc, context=None):
    """Check return code and raise RuntimeError on failure."""
    if rc != QMCKL_SUCCESS:
        msg = _lib.qmckl_string_of_error(rc)
        if msg:
            raise RuntimeError(msg.decode('utf-8'))
        else:
            raise RuntimeError(f"QMCkl error code: {rc}")


# ---------------------------------------------------------------------------
# Pythonic wrapper functions
# ---------------------------------------------------------------------------

# --- Error ---

def string_of_error(error):
    """Convert an exit code into a human-readable string."""
    result = _lib.qmckl_string_of_error(error)
    if result is None:
        return "Unknown error"
    return result.decode('utf-8')


# --- Context ---

def context_create():
    """Create a new QMCkl context.

    Returns:
        int: A new context handle.
    """
    ctx = _lib.qmckl_context_create()
    if ctx == 0:
        raise RuntimeError("Failed to create QMCkl context")
    return ctx


def context_destroy(context):
    """Destroy a QMCkl context.

    Args:
        context: The context to destroy.
    """
    rc = _lib.qmckl_context_destroy(context)
    _check_rc(rc, context)


def context_copy(context):
    """Copy a QMCkl context.

    Args:
        context: The context to copy.

    Returns:
        int: A new context handle.
    """
    new_ctx = _lib.qmckl_context_copy(context)
    if new_ctx == 0:
        raise RuntimeError("Failed to copy QMCkl context")
    return new_ctx


# --- TREXIO ---

def trexio_read(context, file_name):
    """Read a molecular system from a TREXIO file.

    Args:
        context: The QMCkl context.
        file_name (str): Path to the TREXIO file.
    """
    if isinstance(file_name, str):
        file_name = file_name.encode('utf-8')
    rc = _lib.qmckl_trexio_read(context, file_name, len(file_name))
    _check_rc(rc, context)


# --- Numerical precision ---

def set_numprec_precision(context, precision):
    """Set numerical precision."""
    rc = _lib.qmckl_set_numprec_precision(context, precision)
    _check_rc(rc, context)


def set_numprec_range(context, range_val):
    """Set numerical range."""
    rc = _lib.qmckl_set_numprec_range(context, range_val)
    _check_rc(rc, context)


# --- Electron ---

def set_electron_num(context, up_num, down_num):
    """Set the number of up-spin and down-spin electrons."""
    rc = _lib.qmckl_set_electron_num(context, up_num, down_num)
    _check_rc(rc, context)


def set_electron_coord(context, transp, walk_num, coord):
    """Set electron coordinates.

    Args:
        context: The QMCkl context.
        transp (str): Transposition flag ('N' or 'T').
        walk_num (int): Number of walkers.
        coord: Electron coordinates as a list or NumPy array.
    """
    coord = np.asarray(coord, dtype=np.float64)
    if not coord.flags['C_CONTIGUOUS']:
        coord = np.ascontiguousarray(coord)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_set_electron_coord(
        context, transp, walk_num,
        coord.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coord.size
    )
    _check_rc(rc, context)


def get_electron_num(context):
    """Get the total number of electrons."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_electron_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_electron_up_num(context):
    """Get the number of up-spin electrons."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_electron_up_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_electron_down_num(context):
    """Get the number of down-spin electrons."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_electron_down_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_electron_walk_num(context):
    """Get the number of walkers."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_electron_walk_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_electron_coord(context, transp, size_max):
    """Get electron coordinates as a NumPy array.

    Args:
        context: The QMCkl context.
        transp (str): Transposition flag ('N' or 'T').
        size_max (int): Size of the output array.

    Returns:
        numpy.ndarray: Electron coordinates.
    """
    result = np.empty(size_max, dtype=np.float64)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_get_electron_coord(
        context, transp,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


# --- Nucleus ---

def get_nucleus_num(context):
    """Get the number of nuclei."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_nucleus_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_nucleus_charge(context, size_max):
    """Get nuclear charges as a NumPy array."""
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_nucleus_charge(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


def get_nucleus_coord(context, transp, size_max):
    """Get nuclear coordinates as a NumPy array."""
    result = np.empty(size_max, dtype=np.float64)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_get_nucleus_coord(
        context, transp,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


def set_nucleus_num(context, num):
    """Set the number of nuclei."""
    rc = _lib.qmckl_set_nucleus_num(context, num)
    _check_rc(rc, context)


def set_nucleus_charge(context, charge):
    """Set nuclear charges."""
    charge = np.asarray(charge, dtype=np.float64)
    if not charge.flags['C_CONTIGUOUS']:
        charge = np.ascontiguousarray(charge)
    rc = _lib.qmckl_set_nucleus_charge(
        context,
        charge.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        charge.size
    )
    _check_rc(rc, context)


def set_nucleus_coord(context, transp, coord):
    """Set nuclear coordinates."""
    coord = np.asarray(coord, dtype=np.float64)
    if not coord.flags['C_CONTIGUOUS']:
        coord = np.ascontiguousarray(coord)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_set_nucleus_coord(
        context, transp,
        coord.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coord.size
    )
    _check_rc(rc, context)


def get_nucleus_repulsion(context):
    """Get the nuclear repulsion energy."""
    energy = ctypes.c_double()
    rc = _lib.qmckl_get_nucleus_repulsion(context, ctypes.byref(energy))
    _check_rc(rc, context)
    return energy.value


# --- Point ---

def get_point_num(context):
    """Get the number of points."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_point_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def set_point(context, transp, num, coord):
    """Set point coordinates."""
    coord = np.asarray(coord, dtype=np.float64)
    if not coord.flags['C_CONTIGUOUS']:
        coord = np.ascontiguousarray(coord)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_set_point(
        context, transp, num,
        coord.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coord.size
    )
    _check_rc(rc, context)


def get_point(context, transp, size_max):
    """Get point coordinates as a NumPy array."""
    result = np.empty(size_max, dtype=np.float64)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_get_point(
        context, transp,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


# --- AO basis ---

def get_ao_basis_type(context):
    """Get the AO basis type.

    Returns:
        str: Basis type (e.g. 'G' for Gaussian).
    """
    buf = ctypes.create_string_buffer(2)
    rc = _lib.qmckl_get_ao_basis_type(context, buf)
    _check_rc(rc, context)
    return buf.value.decode('utf-8')


def set_ao_basis_type(context, basis_type):
    """Set the AO basis type."""
    if isinstance(basis_type, str):
        basis_type = basis_type.encode('utf-8')
    rc = _lib.qmckl_set_ao_basis_type(context, basis_type, len(basis_type))
    _check_rc(rc, context)


def set_ao_basis_shell_num(context, shell_num):
    """Set the number of shells."""
    rc = _lib.qmckl_set_ao_basis_shell_num(context, shell_num)
    _check_rc(rc, context)


def set_ao_basis_prim_num(context, prim_num):
    """Set the number of primitives."""
    rc = _lib.qmckl_set_ao_basis_prim_num(context, prim_num)
    _check_rc(rc, context)


def set_ao_basis_ao_num(context, ao_num):
    """Set the number of atomic orbitals."""
    rc = _lib.qmckl_set_ao_basis_ao_num(context, ao_num)
    _check_rc(rc, context)


def set_ao_basis_nucleus_index(context, nucleus_index):
    """Set the nucleus index array."""
    nucleus_index = np.asarray(nucleus_index, dtype=np.int64)
    if not nucleus_index.flags['C_CONTIGUOUS']:
        nucleus_index = np.ascontiguousarray(nucleus_index)
    rc = _lib.qmckl_set_ao_basis_nucleus_index(
        context,
        nucleus_index.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
        nucleus_index.size
    )
    _check_rc(rc, context)


def set_ao_basis_nucleus_shell_num(context, nucleus_shell_num):
    """Set the number of shells per nucleus."""
    nucleus_shell_num = np.asarray(nucleus_shell_num, dtype=np.int64)
    if not nucleus_shell_num.flags['C_CONTIGUOUS']:
        nucleus_shell_num = np.ascontiguousarray(nucleus_shell_num)
    rc = _lib.qmckl_set_ao_basis_nucleus_shell_num(
        context,
        nucleus_shell_num.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
        nucleus_shell_num.size
    )
    _check_rc(rc, context)


def set_ao_basis_shell_ang_mom(context, shell_ang_mom):
    """Set the angular momentum of each shell."""
    shell_ang_mom = np.asarray(shell_ang_mom, dtype=np.int32)
    if not shell_ang_mom.flags['C_CONTIGUOUS']:
        shell_ang_mom = np.ascontiguousarray(shell_ang_mom)
    rc = _lib.qmckl_set_ao_basis_shell_ang_mom(
        context,
        shell_ang_mom.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        shell_ang_mom.size
    )
    _check_rc(rc, context)


def set_ao_basis_shell_prim_num(context, shell_prim_num):
    """Set the number of primitives per shell."""
    shell_prim_num = np.asarray(shell_prim_num, dtype=np.int64)
    if not shell_prim_num.flags['C_CONTIGUOUS']:
        shell_prim_num = np.ascontiguousarray(shell_prim_num)
    rc = _lib.qmckl_set_ao_basis_shell_prim_num(
        context,
        shell_prim_num.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
        shell_prim_num.size
    )
    _check_rc(rc, context)


def set_ao_basis_shell_prim_index(context, shell_prim_index):
    """Set the primitive index for each shell."""
    shell_prim_index = np.asarray(shell_prim_index, dtype=np.int64)
    if not shell_prim_index.flags['C_CONTIGUOUS']:
        shell_prim_index = np.ascontiguousarray(shell_prim_index)
    rc = _lib.qmckl_set_ao_basis_shell_prim_index(
        context,
        shell_prim_index.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
        shell_prim_index.size
    )
    _check_rc(rc, context)


def set_ao_basis_shell_factor(context, shell_factor):
    """Set the shell normalization factors."""
    shell_factor = np.asarray(shell_factor, dtype=np.float64)
    if not shell_factor.flags['C_CONTIGUOUS']:
        shell_factor = np.ascontiguousarray(shell_factor)
    rc = _lib.qmckl_set_ao_basis_shell_factor(
        context,
        shell_factor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shell_factor.size
    )
    _check_rc(rc, context)


def set_ao_basis_exponent(context, exponent):
    """Set the primitive exponents."""
    exponent = np.asarray(exponent, dtype=np.float64)
    if not exponent.flags['C_CONTIGUOUS']:
        exponent = np.ascontiguousarray(exponent)
    rc = _lib.qmckl_set_ao_basis_exponent(
        context,
        exponent.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        exponent.size
    )
    _check_rc(rc, context)


def set_ao_basis_coefficient(context, coefficient):
    """Set the primitive coefficients."""
    coefficient = np.asarray(coefficient, dtype=np.float64)
    if not coefficient.flags['C_CONTIGUOUS']:
        coefficient = np.ascontiguousarray(coefficient)
    rc = _lib.qmckl_set_ao_basis_coefficient(
        context,
        coefficient.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coefficient.size
    )
    _check_rc(rc, context)


def set_ao_basis_prim_factor(context, prim_factor):
    """Set the primitive normalization factors."""
    prim_factor = np.asarray(prim_factor, dtype=np.float64)
    if not prim_factor.flags['C_CONTIGUOUS']:
        prim_factor = np.ascontiguousarray(prim_factor)
    rc = _lib.qmckl_set_ao_basis_prim_factor(
        context,
        prim_factor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        prim_factor.size
    )
    _check_rc(rc, context)


def set_ao_basis_ao_factor(context, ao_factor):
    """Set the AO normalization factors."""
    ao_factor = np.asarray(ao_factor, dtype=np.float64)
    if not ao_factor.flags['C_CONTIGUOUS']:
        ao_factor = np.ascontiguousarray(ao_factor)
    rc = _lib.qmckl_set_ao_basis_ao_factor(
        context,
        ao_factor.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ao_factor.size
    )
    _check_rc(rc, context)


# --- MO basis ---

def get_mo_basis_mo_num(context):
    """Get the number of molecular orbitals."""
    mo_num = ctypes.c_int64()
    rc = _lib.qmckl_get_mo_basis_mo_num(context, ctypes.byref(mo_num))
    _check_rc(rc, context)
    return mo_num.value


def set_mo_basis_mo_num(context, mo_num):
    """Set the number of molecular orbitals."""
    rc = _lib.qmckl_set_mo_basis_mo_num(context, mo_num)
    _check_rc(rc, context)


def get_mo_basis_coefficient(context, size_max):
    """Get the MO coefficients as a NumPy array."""
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_mo_basis_coefficient(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


def set_mo_basis_coefficient(context, coefficient):
    """Set the MO coefficients."""
    coefficient = np.asarray(coefficient, dtype=np.float64)
    if not coefficient.flags['C_CONTIGUOUS']:
        coefficient = np.ascontiguousarray(coefficient)
    rc = _lib.qmckl_set_mo_basis_coefficient(
        context,
        coefficient.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coefficient.size
    )
    _check_rc(rc, context)


def get_mo_basis_mo_vgl(context, size_max):
    """Get MO values, gradients and Laplacians as a NumPy array.

    Args:
        context: The QMCkl context.
        size_max (int): Size of the output array (5 * walk_num * elec_num * mo_num).

    Returns:
        numpy.ndarray: MO values, gradients and Laplacians.
    """
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_mo_basis_mo_vgl(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


def get_mo_basis_mo_vgl_inplace(context, size_max):
    """Get MO values, gradients and Laplacians in-place as a NumPy array.

    Args:
        context: The QMCkl context.
        size_max (int): Size of the output array.

    Returns:
        numpy.ndarray: MO values, gradients and Laplacians.
    """
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_mo_basis_mo_vgl_inplace(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        size_max
    )
    _check_rc(rc, context)
    return result


# --- Jastrow CHAMP ---

def set_jastrow_champ_rescale_factor_ee(context, kappa_ee):
    """Set the electron-electron rescale factor for the Jastrow factor."""
    rc = _lib.qmckl_set_jastrow_champ_rescale_factor_ee(context, kappa_ee)
    _check_rc(rc, context)


def set_jastrow_champ_rescale_factor_en(context, kappa_en):
    """Set the electron-nucleus rescale factors for the Jastrow factor."""
    kappa_en = np.asarray(kappa_en, dtype=np.float64)
    if not kappa_en.flags['C_CONTIGUOUS']:
        kappa_en = np.ascontiguousarray(kappa_en)
    rc = _lib.qmckl_set_jastrow_champ_rescale_factor_en(
        context,
        kappa_en.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        kappa_en.size
    )
    _check_rc(rc, context)


def set_jastrow_champ_aord_num(context, aord_num):
    """Set the number of a-order parameters."""
    rc = _lib.qmckl_set_jastrow_champ_aord_num(context, aord_num)
    _check_rc(rc, context)


def set_jastrow_champ_bord_num(context, bord_num):
    """Set the number of b-order parameters."""
    rc = _lib.qmckl_set_jastrow_champ_bord_num(context, bord_num)
    _check_rc(rc, context)


def set_jastrow_champ_cord_num(context, cord_num):
    """Set the number of c-order parameters."""
    rc = _lib.qmckl_set_jastrow_champ_cord_num(context, cord_num)
    _check_rc(rc, context)


def get_jastrow_champ_aord_num(context):
    """Get the number of a-order parameters."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_jastrow_champ_aord_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_jastrow_champ_bord_num(context):
    """Get the number of b-order parameters."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_jastrow_champ_bord_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_jastrow_champ_cord_num(context):
    """Get the number of c-order parameters."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_jastrow_champ_cord_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_jastrow_champ_type_nucl_num(context):
    """Get the number of nucleus types for the Jastrow factor."""
    num = ctypes.c_int64()
    rc = _lib.qmckl_get_jastrow_champ_type_nucl_num(context, ctypes.byref(num))
    _check_rc(rc, context)
    return num.value


def get_jastrow_champ_type_nucl_vector(context, size_max):
    """Get the nucleus type vector."""
    result = np.empty(size_max, dtype=np.int64)
    rc = _lib.qmckl_get_jastrow_champ_type_nucl_vector(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
        size_max
    )
    _check_rc(rc, context)
    return result


# --- Distances / Potentials ---

def get_electron_ee_distance(context, size_max):
    """Get electron-electron distances."""
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_electron_ee_distance(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    _check_rc(rc, context)
    return result


def get_electron_en_distance(context, size_max):
    """Get electron-nucleus distances."""
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_electron_en_distance(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    _check_rc(rc, context)
    return result


def get_electron_ee_potential(context, size_max):
    """Get the electron-electron potential."""
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_electron_ee_potential(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    _check_rc(rc, context)
    return result


def get_electron_en_potential(context, size_max):
    """Get the electron-nucleus potential."""
    result = np.empty(size_max, dtype=np.float64)
    rc = _lib.qmckl_get_electron_en_potential(
        context,
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    )
    _check_rc(rc, context)
    return result
