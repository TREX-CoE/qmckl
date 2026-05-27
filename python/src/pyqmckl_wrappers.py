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
# Hand-written wrapper overrides
#
# These override auto-generated wrappers for functions that need special
# handling (string parameters, context management, etc.).
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


# --- Functions with transposition parameter ---

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


# --- AO basis type (string parameter) ---

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


# --- Single-point update (char parameter) ---

def set_single_point(context, transp, index, coord):
    """Set coordinates of a single point.

    Args:
        context: The QMCkl context.
        transp (str): Transposition flag ('N' or 'T').
        index (int): Index of the point to set.
        coord: Coordinates as a list or NumPy array.
    """
    coord = np.asarray(coord, dtype=np.float64)
    if not coord.flags['C_CONTIGUOUS']:
        coord = np.ascontiguousarray(coord)
    if isinstance(transp, str):
        transp = transp.encode('utf-8')
    rc = _lib.qmckl_set_single_point(
        context, transp, index,
        coord.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coord.size
    )
    _check_rc(rc, context)
