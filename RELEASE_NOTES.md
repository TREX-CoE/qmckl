# QMCkl v1.1.0 Release

We are pleased to announce the release of QMCkl version 1.1.0, the Quantum Monte Carlo kernel library from the TREX Center of Excellence.

## 🎉 Major New Features

### Forces Module for Molecular Dynamics
This release introduces a comprehensive **forces module** that enables:
- **Geometry optimization**: Calculate forces on nuclei (interatomic forces) for finding optimal molecular structures
- **Molecular dynamics simulations**: Time evolution of quantum systems
- **Wavefunction derivatives**: Accurate derivatives with respect to nuclear positions
- **Nonlocal pseudopotential support**: Full force calculation support including nonlocal effects

The forces module implements derivatives of both the Jastrow factor and molecular orbitals, providing the foundation for advanced QMC simulations.

### Single-Electron Move Jastrow Calculations
A new **single-precision Jastrow module** (`qmckl_jastrow_champ_single`) has been added for efficient single-electron updates:
- Computes the difference δJ = J(r', R) - J(r, R) for single-electron moves
- Significantly improves performance for Metropolis-Hastings sampling in QMC
- Includes values, gradients, and Laplacians optimized for single-electron updates
- Essential for efficient Quantum Monte Carlo calculations

## 🔧 Compiler and Build System Enhancements

### Modern Intel Compiler Support
- **Intel ifx**: Support for Intel's next-generation Fortran compiler (oneAPI)
- **Intel icx**: Support for Intel's next-generation C compiler (oneAPI)
- **Improved flags**: Updated from `-Ofast` to `-O3` for better numerical stability
- **MKL updates**: Migrated from `-mkl` to `-qmkl` (modern syntax)
- **Automatic library detection**: Better BLAS/LAPACK handling with Intel compilers

### Enhanced Build Configuration
- **Sanitizer support**: New `--enable-sanitizer` option for comprehensive debugging
  - Address sanitizer
  - Undefined behavior sanitizer
  - Memory leak detection
  - Bounds checking
- **Library flexibility**: Build both shared and static libraries simultaneously
- **Improved portability**: Better compatibility across different platforms and compilers

## 🚀 Continuous Integration & Quality Assurance

### Security and Testing
- **Flawfinder integration**: Automatic security vulnerability scanning
- **Binary distribution workflow**: Automated creation of pre-built packages
- **Updated CI/CD**: Migrated to latest GitHub Actions for better reliability
- **Enhanced testing**: Improved artifact collection for debugging failed builds

### Developer Experience
- **Copilot instructions**: Comprehensive AI-assisted development guidelines
- **Better documentation**: Improved literate programming organization
- **Enhanced examples**: More test data and example calculations

## 📊 What's Changed Since v1.0.0

### New Modules
- `org/qmckl_forces.org` - Complete forces calculation implementation
- `org/qmckl_jastrow_champ_single.org` - Single-electron move Jastrow functions

### Build System
- Support for Intel ifx and icx compilers
- `--enable-sanitizer` configure option
- Both shared and static library support
- Improved CI/CD workflows

### Code Quality
- Security scanning with Flawfinder
- Enhanced test coverage
- Better error handling and debugging support

## 📥 Installation

### From Source Tarball (End Users)
```bash
./configure
make -j
make check
sudo make install
```

### From Repository (Developers)
```bash
./autogen.sh
./configure --prefix=$PWD/_install
make -j
make check
```

### With Intel Compilers
```bash
./configure --with-ifx --with-icx --enable-hpc
make -j
```

## 🔗 Resources

- **Documentation**: https://trex-coe.github.io/qmckl/index.html
- **Source Code**: https://github.com/TREX-CoE/qmckl
- **Issue Tracker**: https://github.com/TREX-CoE/qmckl/issues
- **TREX Project**: https://trex-coe.eu

## 🙏 Acknowledgments

This work is part of the TREX (Targeting Real Chemical Accuracy at the Exascale) project, which has received funding from the European Union's Horizon 2020 Research and Innovation program under grant agreement no. 952165.

## 📝 Full Changelog

**Full Changelog**: https://github.com/TREX-CoE/qmckl/compare/v1.0.0...v1.1.0

---

For detailed information about all changes, see the [NEWS](NEWS) file.
