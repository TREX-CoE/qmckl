
# Python API of the QMCkl library

## Requirements

- `numpy`
- `SWIG` (>= 4.0)

## Manual installation

1. Install the QMCkl as usual
2. Get the latest `qmckl.h` file 
3. `python process.py` to generate `pyqmckl_include.i` list of SWIG patterns
4. `swig -python -py3 -o pyqmckl_wrap.c pyqmckl.i` to generate the SWIG wrapper code in C and `pyqmckl.py` module in Python. 
**Note:** for this to work three files have to be present in the working directory: `pyqmckl.i`, `pyqmckl_include.i` and `numpy.i`.
5. `<c-compiler> -I/usr/include/python3.8 -c -fPIC pyqmckl_wrap.c` to compile the wrapper code into an object file using the `<c-compiler>` (replace with your C compiler, e.g. `gcc`) on your machine
6. `<c-compiler> -shared pyqmckl_wrap.o -lqmckl -o _pyqmckl.so` to produce the final C extension (this requires the `qmckl` library to be installed and present in the linking paths together with all its dependencies like `trexio`)
7. Put the produced `_pyqmckl.so` and `pyqmckl.py` files in the working directory and then run `import pyqmckl`

## Python-ic installation


