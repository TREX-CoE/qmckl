
# Python API of the QMCkl library

This package provides a Python interface to the QMCkl (Quantum Monte Carlo
Kernel Library) using `ctypes`. It is a pure Python package — no compilation
step is needed.


## Requirements

- Python >= 3.6
- `numpy`
- The QMCkl shared library (`libqmckl.so`) installed on the system


## Installation

1. Install the QMCkl C library (see [upstream instructions](https://github.com/TREX-CoE/qmckl))
2. Install the Python package:

```bash
pip install .
```

The package loads `libqmckl.so` at runtime using `ctypes`. Make sure the
library can be found by the dynamic linker. You can either:

- Install QMCkl to a standard system path, or
- Set the `QMCKL_LIBDIR` environment variable to the directory containing
  `libqmckl.so`, or
- Add the directory to `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (macOS).


## Usage

```python
import qmckl as pq

# Create a context
ctx = pq.context_create()

# Read a molecular system from a TREXIO file
pq.trexio_read(ctx, 'molecule.h5')

# Query properties
mo_num = pq.get_mo_basis_mo_num(ctx)
print(f"Number of MOs: {mo_num}")
```


## Running the tests

```bash
cd test/
python test_api.py
```

