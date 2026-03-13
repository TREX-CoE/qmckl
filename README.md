# QMCkl: Quantum Monte Carlo Kernel Library

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

![Build Status](https://github.com/TREX-CoE/qmckl/actions/workflows/test-build.yml/badge.svg?branch=master)

The domain of quantum chemistry needs a library in which the main
kernels of Quantum Monte Carlo (QMC) methods are implemented. In the
library proposed in this project, we expose the main algorithms in a
simple language and provide a standard API and tests to enable the
development of high-performance QMCkl implementations taking
advantage of modern hardware.

See the [source code](https://github.com/TREX-CoE/qmckl/blob/master/org/qmckl.org)
to read the documentation.

## Getting the Source Code

To clone the repository, use:
```bash
git clone https://github.com/TREX-CoE/qmckl.git
```

**Note for Users**: The simplest way to obtain QMCkl is to download a source distribution package.
This repository is primarily for maintainers who write the kernels in org-mode files and 
generate the source code and documentation from these files.

**About Org-mode**: Org-mode is a plain-text format, editable in any text editor—including 
VS Code (which offers excellent support) or Vim (via plugins). Emacs is used here not as a 
text editor but as a command-line tool to generate code from Org-mode files, much like cmake 
or autoconf. Contributors need not master Emacs; they can edit Org-mode files in their 
preferred environment, and our Makefile automates the rest.

# Prerequisites and Dependencies

## Required Dependencies

Before building QMCkl, you should install the following dependencies:

### 1. TREXIO Library (Highly Recommended)

**TREXIO** (TREX I/O) is a library for reading and writing wave function data files.
While TREXIO is enabled by default and highly encouraged, it can be disabled during 
configuration with `./configure --without-trexio` if needed.

- **Repository**: https://github.com/TREX-CoE/trexio
- **Releases**: https://github.com/TREX-CoE/trexio/releases
- **Installation** (from latest release):
  ```bash
  # Download the latest tar.gz from the releases page, then:
  tar -xzf trexio-*.tar.gz
  cd trexio-*
  ./configure --prefix=/path/to/install
  make
  make check
  sudo make install
  ```

- **Setting PKG_CONFIG_PATH**: If you installed TREXIO in a non-standard location 
  (i.e., not `/usr` or `/usr/local`), you need to set the `PKG_CONFIG_PATH` 
  environment variable:
  ```bash
  export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH
  ```

  Alternatively, you can specify the TREXIO location during configuration:
  ```bash
  ./configure --with-trexio=/path/to/install
  ```

  Or set the environment variables directly:
  ```bash
  export TREXIO_CFLAGS="-I/path/to/install/include"
  export TREXIO_LIBS="-L/path/to/install/lib -ltrexio"
  ```

### 2. BLAS and LAPACK

**BLAS** (Basic Linear Algebra Subprograms) and **LAPACK** (Linear Algebra PACKage) 
are required for linear algebra operations.

**Performance Note**: Generic BLAS implementations will result in low performance. 
Vendor-specific BLAS libraries such as Intel MKL or ARM Performance Libraries (ARMPL) 
are highly recommended for optimal performance.

- **Installation on Ubuntu/Debian**:
  ```bash
  sudo apt-get install libblas-dev liblapack-dev
  ```

- **Installation on RHEL/CentOS/Fedora**:
  ```bash
  sudo yum install blas-devel lapack-devel
  ```

- **Recommended**: Intel MKL (Math Kernel Library) or ARM Performance Libraries provide 
  optimized BLAS/LAPACK. When using Intel compilers with `--with-icc` or `--with-ifort`, 
  MKL is automatically used.

### 3. HDF5 (Optional, Highly Recommended)

**HDF5** is an optional dependency of TREXIO for data storage, but it is highly recommended 
and enabled by default in TREXIO's configuration.

- **Installation on Ubuntu/Debian**:
  ```bash
  sudo apt-get install libhdf5-dev
  ```

- **Installation on RHEL/CentOS/Fedora**:
  ```bash
  sudo yum install hdf5-devel
  ```

### 4. Build Tools

- **Autotools**: `autoconf`, `automake`, `libtool`
- **C Compiler**: GCC or Intel ICC
- **Fortran Compiler**: gfortran or Intel ifort/ifx
- **pkg-config** (optional): Helps find library dependencies

- **Installation on Ubuntu/Debian**:
  ```bash
  sudo apt-get install build-essential autoconf automake libtool pkg-config gfortran
  ```

- **Installation on RHEL/CentOS/Fedora**:
  ```bash
  sudo yum install gcc gcc-gfortran autoconf automake libtool pkgconfig
  ```

## Optional Dependencies

### For Maintainers (Developer Mode)

If you are working with the org-mode source files (maintainer mode):

- **Emacs** (>= 26): For tangling org-mode files
  ```bash
  sudo apt-get install emacs  # Ubuntu/Debian
  sudo yum install emacs      # RHEL/CentOS/Fedora
  ```

- **Python 3**: For build scripts
- **AWK**: Usually pre-installed

### For Python API

- **SWIG** (>= 4.0, optional): For generating Python bindings
  ```bash
  sudo apt-get install swig   # Ubuntu/Debian
  sudo yum install swig       # RHEL/CentOS/Fedora
  ```

# Installation

## Quick Start for Users

**Prerequisites**: Ensure all [required dependencies](#required-dependencies) are installed first.

1. **Download a source distribution**:
   - Get the latest release from: https://github.com/TREX-CoE/qmckl/releases
   ```bash
   tar -xzf qmckl-*.tar.gz
   cd qmckl-*
   ```

2. **Configure the build**:
   ```bash
   ./configure
   ```

   If TREXIO or other dependencies are in non-standard locations:
   ```bash
   export PKG_CONFIG_PATH=/path/to/trexio/lib/pkgconfig:$PKG_CONFIG_PATH
   ./configure
   ```

3. **Build and test**:
   ```bash
   make -j
   make -j check
   ```

4. **Install**:
   ```bash
   sudo make install
   sudo make installcheck
   ```

## Advanced Build Options

### Out-of-Source Builds (Recommended)

Out-of-source builds keep your source directory clean and allow multiple 
configurations simultaneously:

```bash
mkdir _build
cd _build
../configure --prefix=$PWD/_install
make -j
make -j check
make install
```

### Debug Build

For development and debugging:

```bash
./configure --enable-debug --prefix=/path/to/install
make -j
make -j check
```

### Optimized Build with GCC

For high-performance builds with GCC:

```bash
./configure \
  CC=gcc \
  CFLAGS="-g -O2 -march=native -flto -fno-trapping-math -fno-math-errno -ftree-vectorize" \
  FC=gfortran \
  FCFLAGS="-g -O2 -march=native -flto -ftree-vectorize" \
  --enable-hpc
make -j
make -j check
```

### Optimized Build with Intel Compilers

For high-performance builds with Intel compilers (which include MKL for BLAS/LAPACK):

```bash
./configure \
  --with-icc \
  --with-ifort \
  --enable-hpc \
  --prefix=/path/to/install
make -j
make -j check
```

Alternatively, with newer Intel compilers:

```bash
./configure \
  --with-icx \
  --with-ifx \
  --enable-hpc \
  --prefix=/path/to/install
make -j
make -j check
```

## For Maintainers (Developer Mode)

If you are working with the org-mode source files (when `.maintainer_mode` file exists):

**Additional Prerequisites**: Emacs (>= 26), Python 3, AWK

```bash
./autogen.sh
./configure --prefix=$PWD/_install
make
make check
```

### Multiple Build Configurations

You can create different build configurations simultaneously:

```bash
./autogen.sh

# Debug build with GCC
mkdir -p _build_gcc_debug/_install
cd _build_gcc_debug
../configure --enable-debug --prefix=$PWD/_install
make -j
make -j check
cd ..

# Optimized build with Intel compilers
mkdir -p _build_intel_fast/_install
cd _build_intel_fast
../configure --prefix=$PWD/_install --enable-hpc --with-icc --with-ifort 
make -j
make -j check
cd ..
```

## Python API

The Python API allows you to use QMCkl from Python.

### Prerequisites

- **SWIG** (>= 4.0) is required
- QMCkl shared library must be installed first

### Installation

1. **Install the QMCkl C library** following the instructions above

2. **Install the Python package**:
   ```bash
   make python-install
   ```

3. **Test the installation**:
   ```bash
   make python-test
   ```

### Usage Example

See [test_api.py](https://github.com/TREX-CoE/qmckl/blob/master/python/test/test_api.py) 
for examples of using the Python API.

**Recommendation**: Use [virtual environments](https://docs.python.org/3/tutorial/venv.html) 
to avoid compatibility issues:

```bash
python3 -m venv qmckl_env
source qmckl_env/bin/activate
make python-install
make python-test
```

## Linking to Your Program

After installation, link your program against QMCkl by adding `-lqmckl` to your compiler options.

### Setting Library Paths

If QMCkl was installed with a custom prefix, update your library and include search paths:

```bash
export LIBRARY_PATH=$LIBRARY_PATH:/path/to/qmckl/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/qmckl/lib
export CPATH=$CPATH:/path/to/qmckl/include
```

### Using CMake

If your project uses CMake, you can use the provided 
[FindQMCKL.cmake](https://github.com/TREX-CoE/qmckl/blob/master/cmake/FindQMCKL.cmake) 
module to automatically find and link QMCkl:

```cmake
list(APPEND CMAKE_MODULE_PATH "/path/to/qmckl/cmake")
find_package(QMCKL REQUIRED)
target_link_libraries(your_target QMCKL::qmckl)
```

## Installation with GNU Guix

QMCkl can be installed using the [GNU Guix](https://guix.gnu.org) package manager.

The [qmckl.scm](https://github.com/TREX-CoE/qmckl/blob/master/tools/qmckl.scm) 
file contains the package definition:

```bash
guix package \
  --profile=$GUIX_PROFILE \
  --load-path=<path_to_trexio_scm> \
  --cores=<n_cores> \
  --install-from-file=qmckl.scm
```

Where `<path_to_trexio_scm>` should point to the directory containing TREXIO's 
[trexio.scm](https://github.com/TREX-CoE/trexio/blob/master/tools/trexio.scm) 
file (e.g., `~/trexio/tools/`).

You can choose between development (`qmckl-dev`) and stable (`qmckl-hpc`) versions 
by modifying the return value in `qmckl.scm`.

# Troubleshooting

## Common Issues

### TREXIO not found

**Error message**:
```
configure: error: Package requirements (trexio) were not met:
Package 'trexio', required by 'virtual:world', not found
```

**Solution**: 
1. Ensure TREXIO is installed (see [TREXIO installation](#1-trexio-library))
2. If installed in a non-standard location, set `PKG_CONFIG_PATH`:
   ```bash
   export PKG_CONFIG_PATH=/path/to/trexio/lib/pkgconfig:$PKG_CONFIG_PATH
   ```
3. Or use `--with-trexio`:
   ```bash
   ./configure --with-trexio=/path/to/trexio
   ```

### BLAS/LAPACK not found

**Error message**:
```
configure: error: BLAS was not found.
configure: error: LAPACK was not found.
```

**Solution**: Install BLAS and LAPACK development packages:
- **Ubuntu/Debian**: `sudo apt-get install libblas-dev liblapack-dev`
- **RHEL/CentOS**: `sudo yum install blas-devel lapack-devel`

### Emacs required error (Maintainer mode only)

**Error message**:
```
Error: Emacs is required for org-mode.
```

**Solution**: 
- If you're a **user** (not a maintainer), download a source distribution instead of cloning the repository
- If you're a **maintainer**, install Emacs: `sudo apt-get install emacs` or `sudo yum install emacs`

### HDF5 not found (TREXIO dependency)

**Solution**: Install HDF5 development package:
- **Ubuntu/Debian**: `sudo apt-get install libhdf5-dev`
- **RHEL/CentOS**: `sudo yum install hdf5-devel`

## Getting Help

- **Documentation**: https://trex-coe.github.io/qmckl/index.html
- **Issues**: https://github.com/TREX-CoE/qmckl/issues
- **TREX Project**: https://trex-coe.eu

# Advanced Topics

## Verificarlo Support

Verificarlo is a tool for numerical verification. QMCkl can be built with Verificarlo support,
but it is **not a required dependency**.

All Verificarlo functions are only called when explicitly enabled and are ignored otherwise.

To enable Verificarlo support:

```bash
./configure \
  CC="verificarlo-f" \
  FC="verificarlo-f" \
  --prefix=$PWD/_install \
  --enable-vfc_ci \
  --host=x86_64
make -j
make -j check
```

Building without the `--enable-vfc_ci` flag will ignore all Verificarlo-related code.

## Configure Options

Run `./configure --help` to see all available options. Common options include:

- `--prefix=/path/to/install` - Installation directory
- `--enable-hpc` - Enable HPC-optimized functions
- `--enable-debug` - Build with debugging symbols and checks
- `--with-icc` - Use Intel icc C compiler
- `--with-ifort` - Use Intel ifort Fortran compiler
- `--with-icx` - Use Intel icx C compiler (newer)
- `--with-ifx` - Use Intel ifx Fortran compiler (newer)
- `--with-trexio=/path` - Specify TREXIO installation directory
- `--without-trexio` - Build without TREXIO support (not recommended)
- `--with-qmckldgemm=/path` - Use custom QMCKL DGEMM library
- `--without-openmp` - Disable OpenMP support
- `--enable-python` - Enable Python API build
- `--enable-vfc_ci` - Enable Verificarlo support
- `--enable-malloc-trace` - Enable malloc/free debugging
- `--enable-prof` - Enable profiling
- `--with-efence` - Use ElectricFence library for memory debugging

# Supported Platforms

QMCkl has been tested on:

- **CPUs**: x86 and ARM CPUs
- **Operating Systems**: Ubuntu, Debian, RHEL, CentOS, Fedora, macOS (with Homebrew or MacPorts)
- **Compilers**: GCC, Intel ICC/ICX, Intel ifort/ifx

# License

QMCkl is licensed under the **BSD 3-Clause License**. 
See the [LICENSE](LICENSE) file for details.

**Copyright** (c) 2020, TREX Center of Excellence

# Citation

If you use QMCkl in your research, please cite:

[arXiv:2512.16677v1](https://arxiv.org/abs/2512.16677v1)

# Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes (in org-mode files if you're modifying code)
4. Run tests with `make check`
5. Submit a pull request

For major changes, please open an issue first to discuss the proposed changes.

# Acknowledgments

------------------------------

![European flag](https://trex-coe.eu/sites/default/files/inline-images/euflag.jpg)

[TREX: Targeting Real Chemical Accuracy at the Exascale](https://trex-coe.eu) project has received funding from the European Union's Horizon 2020 - Research and Innovation program - under grant agreement no. 952165. The content of this document does not represent the opinion of the European Union, and the European Union is not responsible for any use that might be made of such content.
