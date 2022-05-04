#!/bin/bash

set -x
set -e

./build_pyqmckl.sh

# copy swig-produced pyqmckl.py module into the pyqmckl/ folder
#cp src/qmckl.py qmckl/
cp src/qmckl.py ./

# install using pip
pip install .

