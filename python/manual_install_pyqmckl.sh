#!/bin/bash

set -x
set -e

./build_pyqmckl.sh

cd src/

# compile the wrapper code
cc -c -fPIC `pkg-config --cflags qmckl` -I/usr/include/python3.8 pyqmckl_wrap.c -o pyqmckl_wrap.o

# link against the previously installed QMCkl library (as detected by pkg-config)
cc -shared pyqmckl_wrap.o `pkg-config --libs qmckl` -o _pyqmckl.so

cd ..

# test
cp src/_pyqmckl.so src/pyqmckl.py test/

cd test/ 
python test_api.py
cd ..
