#!/bin/bash

set -e
set -x

cp ../include/qmckl.h .

# check if qmckl header exists
if [[ ! -f 'qmckl.h' ]]; then
  echo "qmckl.h NOT FOUND"
  exit 1
fi

# process the qmckl header file to get patterns for SWIG
python process.py

# check if SWIG files exist
SWIG_LIST='pyqmckl.i pyqmckl_include.i numpy.i'
for file in $SWIG_LIST; do
  if [[ ! -f $file ]]; then
    echo "$file NOT FOUND"
    exit 1
  fi
done

# run SWIG interface file to produce the Python wrappers
swig -python -py3 -o pyqmckl_wrap.c pyqmckl.i 

# compile the wrapper code
cc -c -fPIC `pkg-config --cflags qmckl` -I/usr/include/python3.8 pyqmckl_wrap.c -o pyqmckl_wrap.o

# link against the previously installed QMCkl library (as detected by pkg-config)
cc -shared pyqmckl_wrap.o `pkg-config --libs qmckl` -o _pyqmckl.so

# test
cp _pyqmckl.so pyqmckl.py -- test/
cd test 
python test_api.py

