#!/bin/bash

python ./tools/build_makefile.py
autoreconf -i -Wall --no-recursive
