#!/bin/bash

export srcdir="."
python ${srcdir}/tools/build_makefile.py
autoreconf -i -Wall --no-recursive
