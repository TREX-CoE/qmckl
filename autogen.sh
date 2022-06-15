#!/bin/bash

export srcdir="."
python3 ${srcdir}/tools/build_makefile.py
autoreconf -i -Wall --no-recursive
