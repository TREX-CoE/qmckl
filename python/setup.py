#!/usr/bin/env python3
"""
setup.py file for pyqmckl package
"""

from setuptools import setup, Extension
from os.path import join


# Read the long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# Read the version string from the file
VERSIONFILE = "pyqmckl/_version.py"
try:
    exec(open(VERSIONFILE).read())
except:
    raise IOError("Could not open the version file %s." % (VERSIONFILE, ))

version_r = __version__
if not version_r:
    raise RuntimeError("Unable to find a version string in %s." % (VERSIONFILE, ))


# Define the name of the Python package
mod_name = 'pyqmckl'


# Define pyqmckl extension module based on SWIG interface file (requires qmckl.h)
pyqmckl_module =  Extension(name               = mod_name + '._' + mod_name,
                            sources            = [ join('src', mod_name + '_wrap.c') ],
                            #include_dirs       = [numpy_includedir],
                            #library_dirs       = [],
                            #runtime_library_dirs = [],
                            libraries          = ['qmckl'],
                            extra_compile_args = ['-Wall'],
                            #extra_link_args    = [h5_ldflags],
                            #swig_opts          = ['-py3' , '-builtin'],
                            depends            = [ join('src', 'qmckl.h') ],
                            language           = 'c'
                            )


setup(name             = mod_name,
      version          = version_r,
      author           = "TREX-CoE",
      author_email     = "posenitskiy@irsamc.ups-tlse.fr",
      description      = """Python API of the QMCkl library""",
      long_description = long_description,
      long_description_content_type = "text/markdown",
      ext_modules      = [pyqmckl_module],
      py_modules       = [mod_name],
      packages         = [mod_name],
      url              = 'https://github.com/TREX-CoE/qmckl',
      license          = 'BSD',
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
      install_requires = ['numpy>=1.17.3']
      )
