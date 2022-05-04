#!/usr/bin/env python3
"""
setup.py file for qmckl package
"""

from setuptools import setup, Extension
from os.path import join


# Read the long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# Define the name of the Python package
MODULE_NAME = "qmckl"

# Define qmckl extension module based on SWIG interface file (requires qmckl.h)
qmckl_module =  Extension(name               = "._" + MODULE_NAME,
                            sources            = [ join("src", MODULE_NAME + "_wrap.c") ],
                            #include_dirs       = [numpy_includedir],
                            #library_dirs       = [],
                            libraries          = ["qmckl"],
                            extra_compile_args = ["-Wall"],
                            #extra_link_args    = [h5_ldflags],
                            depends            = [ join("src", "qmckl.h") ],
                            language           = "c"
                            )


setup(name             = MODULE_NAME,
      version          = "0.2.0",
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
      install_requires = ["numpy>=1.17.3"]
      )
