
# Python API of the QMCkl library


## Requirements

- `setuptools`
- `numpy`
- `swig` (>= 4.0)


## Manual installation

1. Install the QMCkl library (see upstream instructions)
2. `./manual_install_qmckl.sh` which should do the following
3. Copy the produced `_qmckl.so` and `qmckl.py` files into your working directory and do not forget to `import qmckl` in your Python scripts

The second step executes the following under the hood:

1. `./build_qmckl.sh`
2. `<c-compiler> -I/usr/include/python3.8 -c -fPIC qmckl_wrap.c` to compile the wrapper code into an object file using the `<c-compiler>` (replace with your C compiler, e.g. `gcc`) on your machine
3. `<c-compiler> -shared qmckl_wrap.o -lqmckl -o _qmckl.so` to produce the final C extension (this requires the `qmckl` library to be installed and present in the linking paths together with all its dependencies like `trexio`)


## Python-ic installation (recommended)

1. Install the QMCkl library (see upstream instructions)
2. `./pip_install_qmckl.sh`

The last step runs `./build_qmckl.sh`, copies the result into the `qmckl/` directory and 
then runs `pip install .` to install the `qmckl` Python package in your environment.


## SWIG pre-processing

Both aforementioned steps call `build_qmckl.sh` script which does the following pre-processing for SWIG

1. Copy the latest `qmckl.h` file fron `include/` into the `src/` directory
2. `python process_header.py` to generate `qmckl_include.i` list of SWIG patterns
3. `swig -python -py3 -builtin -threads -o qmckl_wrap.c qmckl.i` to generate the SWIG wrapper code in C and `qmckl.py` module in Python. 
**Note:** for this to work three files have to be present in the working directory: `qmckl.i`, `qmckl_include.i` and `numpy.i`.

