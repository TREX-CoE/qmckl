#!/bin/bash

set -x
set -e

# swig pre-processing
./build_qmckl.sh

cd src/

# compile the wrapper code
cc -c -fPIC `pkg-config --cflags qmckl` -I/usr/include/python3.8 qmckl_wrap.c -o qmckl_wrap.o

# link against the previously installed QMCkl library (as detected by pkg-config)
cc -shared qmckl_wrap.o `pkg-config --libs qmckl` -o _qmckl.so

cd ..

# copy the produced files into the test dir
cp src/_qmckl.so src/qmckl.py test/

# run tests
cd test/ 
python test_api.py
cd ..
