#!/bin/sh

export srcdir="."
python3 ${srcdir}/tools/build_makefile.py
autoreconf -vi -Wall --no-recursive
