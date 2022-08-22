#!/usr/bin/env python3
"""
setup.py file for qmckl package
"""

import os, sys
from setuptools import setup, Extension
from os.path import join


# Read the long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# this was recommended to solve the problem of the missing numpy header files
try:
    import numpy
except ImportError:
    raise Exception("numpy Python package cannot be imported.")

numpy_includedir = numpy.get_include()

# Define the name of the Python package
MODULE_NAME = "qmckl"

# derive the QMCkl libdir and includedir
QMCKL_LIBDIR     = os.environ.get("QMCKL_LIBDIR", None)
QMCKL_INCLUDEDIR = os.environ.get("QMCKL_INCLUDEDIR", None)

libdir_undefined     = QMCKL_LIBDIR is None or QMCKL_LIBDIR==""
includedir_undefined = QMCKL_INCLUDEDIR is None or QMCKL_INCLUDEDIR==""


# Define qmckl extension module based on SWIG interface file (requires qmckl.h)
qmckl_module   =  Extension(name                 = "._" + MODULE_NAME,
                            sources              = [ join("src", MODULE_NAME + "_wrap.c") ],
                            include_dirs         = [numpy_includedir, QMCKL_INCLUDEDIR],
                            #library_dirs         = [QMCKL_LIBDIR],
                            runtime_library_dirs = [QMCKL_LIBDIR],
                            libraries            = ["qmckl"],
                            extra_compile_args   = ["-Wall"],
                            extra_link_args      = ["-L" + QMCKL_LIBDIR],
                            depends              = [ join("src", "qmckl.h") ],
                            language             = "c"
                            )


setup(name             = MODULE_NAME,
      version          = "0.3.1",
      author           = "TREX-CoE",
      author_email     = "posenitskiy@irsamc.ups-tlse.fr",
      description      = """Python API of the QMCkl library""",
      long_description = long_description,
      long_description_content_type = "text/markdown",
      ext_modules      = [qmckl_module],
      py_modules       = [MODULE_NAME],
      url              = "https://github.com/TREX-CoE/qmckl",
      license          = "BSD",
      classifiers=[
         "Intended Audience :: Science/Research",
         "Intended Audience :: Developers",
         "Topic :: Scientific/Engineering",
         "Programming Language :: C",
         "Programming Language :: Python",
         "Programming Language :: Python :: 3",
         "Programming Language :: Python :: 3 :: Only",
         "Programming Language :: Python :: Implementation :: CPython",
         "License :: OSI Approved :: BSD License",
         "Operating System :: POSIX",
         "Operating System :: Unix",
         "Operating System :: MacOS"
      ],
      python_requires = ">=3.0",
      # The ABI incompatibility of numpy is a pain, for now set the 
      # min numpy version such that we have wheels for CPython 3.5 & 3.6
      install_requires = ["numpy>=1.13.3"]
      )
