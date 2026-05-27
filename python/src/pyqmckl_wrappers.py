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
