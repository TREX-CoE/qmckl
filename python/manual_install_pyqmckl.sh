#!/bin/bash

set -x
set -e

# swig pre-processing
./build_pyqmckl.sh

cd src/

# compile the wrapper code
cc -c -fPIC `pkg-config --cflags qmckl` -I/usr/include/python3.8 pyqmckl_wrap.c -o pyqmckl_wrap.o

# link against the previously installed QMCkl library (as detected by pkg-config)
cc -shared pyqmckl_wrap.o `pkg-config --libs qmckl` -o _pyqmckl.so

cd ..

# copy the produced files into the test dir
cp src/_pyqmckl.so src/pyqmckl.py test/

# run tests
cd test/ 
python test_api.py
cd ..
