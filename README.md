# QMCkl: Quantum Monte Carlo Kernel Library

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

![Build Status](https://github.com/TREX-CoE/qmckl/workflows/test-build/badge.svg?branch=master)

The domain of quantum chemistry needs a library in which the main
kernels of Quantum Monte Carlo (QMC) methods are implemented. In the
library proposed in this project, we expose the main algorithms in a
simple language and provide a standard API and tests to enable the
development of high-performance QMCkl implementations taking
advantage of modern hardware.

See the [source code](https://github.com/TREX-CoE/qmckl/blob/master/org/qmckl.org)
to read the documentation.


To clone the repository, use:
```
git clone https://github.com/TREX-CoE/qmckl.git
```

# Installation

The simplest way to obtain the source files of QMCkl is to download a source
distribution. This particular repository is for maintainers, who write the kernels
in org-mode files and produce the source code and the documentation from these files.

## For maintainers

```
./autogen.sh
./configure --prefix=$PWD/_install

make
make check
```

## For users

Obtain a source distribution.

To build the documentation version:

```
./configure
```

To build an optimized version with Intel compilers:
```
./configure \
   --with-icc \
   --with-ifort \
   --enable-hpc \
   --with-openmp
```

To build an optimized version with GCC:
```
./configure \
  CC=gcc \
  CFLAGS="-g -O2 -march=native  -flto -fno-trapping-math -fno-math-errno -ftree-vectorize" \
  FC=gfortran \
  FCFLAGS="-g -O2 -march=native  -flto -ftree-vectorize" \
  --enable-hpc \
  --with-openmp
```


Then, compile with:
```
make -j
make -j check
sudo make install
sudo make installcheck
```


## Linking to your program

The `make install` command takes care of installing the QMCkl shared library on the user machine.
Once installed, add `-lqmckl` to the list of compiler options.

In some cases (e.g. when using custom `prefix` during configuration), the QMCkl library might end up installed in a directory, which is absent in the default `$LIBRARY_PATH`.
In order to link the program against QMCkl, the search paths can be modified as follows:

`export LIBRARY_PATH=$LIBRARY_PATH:<path_to_qmckl>/lib`

(same holds for `$LD_LIBRARY_PATH`). The `<path_to_qmckl>` has to be replaced with the prefix used during the installation.

If your project relies on the CMake build system, feel free to use the
[FindQMCKL.cmake](https://github.com/TREX-CoE/qmckl/blob/master/cmake/FindQMCKL.cmake)
module to find and link the QMCkl library automatically.


## Verificarlo CI

Since Verificarlo should not be a dependency of QMCkl, all Verificarlo
functions are called only when the support is explicitely enabled (and ignored
by the preprocessor otherwise). To enable vfc_ci support, the library should be
configured with the following command :

```
./configure \
  CC="verificarlo-f" \
  FC="verificarlo-f" \
  --prefix=$PWD/_install \
  --enable-vfc_ci \
  --host=x86_64 \
```

where CC and FC are set to verificarlo-f, and support is explicitely enabled
with the --enable-vfc_ci flag. Configuring the library with the "standard"
command will cause all calls to Verificarlo related functions to be ignored,
and the library will be built as usual.


------------------------------

![European flag](https://trex-coe.eu/sites/default/files/inline-images/euflag.jpg)
[TREX: Targeting Real Chemical Accuracy at the Exascale](https://trex-coe.eu) project has received funding from the European Union’s Horizon 2020 - Research and Innovation program - under grant agreement no. 952165. The content of this document does not represent the opinion of the European Union, and the European Union is not responsible for any use that might be made of such content.

