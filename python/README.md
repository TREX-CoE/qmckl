
# Python API of the QMCkl library


## Requirements

- `setuptools`
- `numpy`
- `swig` (>= 4.0)


## Manual installation

1. Install the QMCkl library (see upstream instructions)
2. `./manual_install_pyqmckl.sh` which should do the following
3. Copy the produced `_pyqmckl.so` and `pyqmckl.py` files into your working directory and do not forget to `import pyqmckl` in your Python scripts

The second step executes the following under the hood:

1. `./build_pyqmckl.sh`
2. `<c-compiler> -I/usr/include/python3.8 -c -fPIC pyqmckl_wrap.c` to compile the wrapper code into an object file using the `<c-compiler>` (replace with your C compiler, e.g. `gcc`) on your machine
3. `<c-compiler> -shared pyqmckl_wrap.o -lqmckl -o _pyqmckl.so` to produce the final C extension (this requires the `qmckl` library to be installed and present in the linking paths together with all its dependencies like `trexio`)


## Python-ic installation (recommended)

1. Install the QMCkl library (see upstream instructions)
2. `./pip_install_pyqmckl.sh`

The last step runs `./build_pyqmckl.sh`, copies the result into the `pyqmckl/` directory and 
then runs `pip install .` to install the `pyqmckl` Python package in your environment.


## SWIG pre-processing

Both aforementioned steps call `build_pyqmckl.sh` script which does the following pre-processing for SWIG

1. Copy the latest `qmckl.h` file fron `include/` into the `src/` directory
2. `python process_header.py` to generate `pyqmckl_include.i` list of SWIG patterns
3. `swig -python -py3 -builtin -threads -o pyqmckl_wrap.c pyqmckl.i` to generate the SWIG wrapper code in C and `pyqmckl.py` module in Python. 
**Note:** for this to work three files have to be present in the working directory: `pyqmckl.i`, `pyqmckl_include.i` and `numpy.i`.

