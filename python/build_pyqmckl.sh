#!/bin/bash

set -e
set -x

cp ../include/qmckl.h src/

cd src/

# check if qmckl header exists
if [[ ! -f 'qmckl.h' ]]; then
  echo "qmckl.h NOT FOUND"
  exit 1
fi

# process the qmckl header file to get patterns for SWIG
python process_header.py

# check if SWIG files exist
SWIG_LIST='pyqmckl.i pyqmckl_include.i numpy.i'
for file in $SWIG_LIST; do
  if [[ ! -f $file ]]; then
    echo "$file NOT FOUND"
    exit 1
  fi
done

# run SWIG interface file to produce the Python wrappers
# `-threads` to release GIL before the call to C functions (for multi-threading)
# `-builtin` to activate Python-ic built-in types
swig -python -py3 -builtin -o qmckl_wrap.c pyqmckl.i 

cd ..
