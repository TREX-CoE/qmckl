# Copilot Instructions for QMCkl

## Repository Overview

QMCkl (Quantum Monte Carlo kernel library) is a C/Fortran library that implements the main kernels of Quantum Monte Carlo (QMC) methods. The library provides a standard API with implementations optimized for modern hardware.

- **License**: BSD 3-Clause License
- **Copyright**: TREX Center of Excellence
- **Project**: Part of the TREX (Targeting Real Chemical Accuracy at the Exascale) project funded by the European Union's Horizon 2020

## Development Approach: Literate Programming

**IMPORTANT**: This project uses literate programming with Emacs org-mode files.

### Source Code Structure

- **Primary sources**: All code lives in `.org` files in the `/org` directory
- **Generated files**: C and Fortran source files are generated from org-mode files
- **Never edit generated files directly**: Always edit the `.org` files
- **Documentation**: The org-mode files are both documentation AND source code

### Org-Mode File Structure

Each `.org` file contains:
- Detailed documentation with mathematical formulas and explanations
- Code blocks marked with `#+begin_src c`, `#+begin_src f90`, etc.
- The code blocks are tangled (extracted) to produce actual source files
- Headers specify tangling destinations using `:tangle` properties

### Key Org Files

- `org/qmckl.org` - Main documentation and introduction
- `org/qmckl_context.org` - Context management (fundamental to the API)
- `org/qmckl_error.org` - Error handling conventions
- `org/qmckl_memory.org` - Memory management
- `org/qmckl_nucleus.org`, `org/qmckl_electron.org`, etc. - Domain-specific modules

## Build System

### For Maintainers (working with org-mode sources)

```bash
./autogen.sh
./configure --prefix=$PWD/_install
make
make check
```

### Build Requirements for Maintainers

- Emacs >= 26 (for tangling org-mode files)
- Autotools (autoconf, automake)
- Python 3
- C compiler (GCC or Intel ICC)
- Fortran compiler (gfortran or Intel ifort)
- BLAS and LAPACK libraries
- HDF5 library
- TREXIO library (https://github.com/TREX-CoE/trexio)

### Build Modes

- **Debug mode**: `./configure --enable-debug`
- **HPC mode**: `./configure --enable-hpc` (optimized for performance)
- **With Verificarlo**: `./configure --enable-vfc_ci` (for numerical verification)

### Out-of-Source Builds (Recommended)

```bash
mkdir _build
cd _build
../configure [options]
make -j
make -j check
```

## Testing

- Run tests with: `make check`
- Tests are defined within the org-mode files
- Test suite logs are available in `test-suite.log` after running tests
- Python tests: `make python-test` (after installing Python API)

## Coding Conventions

### Error Handling

- **All functions return `qmckl_exit_code`** (an `int32_t` type)
- **Success**: Return `QMCKL_SUCCESS` (value 0)
- **Errors**: Return specific error codes like `QMCKL_INVALID_ARG_1`, `QMCKL_INVALID_ARG_2`, etc.
- **Never abort**: Library functions should never cause the program to abort
- **No I/O**: Library functions should not perform input/output operations directly

### Naming Conventions

- Functions: `qmckl_module_function_name` (e.g., `qmckl_context_create`)
- Types: `qmckl_type_name` (e.g., `qmckl_exit_code`, `qmckl_context`)
- Macros: `QMCKL_MACRO_NAME` (e.g., `QMCKL_SUCCESS`)

### Code Organization

- Each module has its own org-mode file
- Private headers: `qmckl_module_private_type.h`, `qmckl_module_private_func.h`
- Public headers: Generated and installed in `include/`
- Fortran interface: `qmckl_f.f90` provides Fortran bindings

### Language-Specific Guidelines

**C Code:**
- Use standard C99 or later
- Include guards in all headers
- Use `stdint.h` types (e.g., `int32_t`, `int64_t`)
- Thread-safe: Use pthread for thread safety when needed

**Fortran Code:**
- Use Fortran 2008 or later
- Use `iso_c_binding` for C interoperability
- Follow the module structure defined in org files

## Python API

- Built with SWIG (>= 4.0)
- Install with: `make python-install`
- Test with: `make python-test`
- Example usage in `python/test/test_api.py`

## Continuous Integration

The project uses GitHub Actions with multiple build configurations:
- **standard**: Standard build with Python API tests
- **debug**: Debug mode build with additional checks
- **hpc**: High-performance optimized build
- **macos**: macOS compatibility testing

CI workflow is defined in `.github/workflows/test-build.yml`

## Dependencies

### Required
- TREXIO library for reading/writing data files
- BLAS and LAPACK for linear algebra
- HDF5 for data storage

### Optional
- SWIG for Python bindings
- Verificarlo for numerical verification

## Key Files and Directories

- `/org` - Literate programming source files (org-mode)
- `/src` - Generated C source files (do not edit directly)
- `/include` - Public headers
- `/python` - Python API bindings
- `/tools` - Build tools and utilities
- `/share` - Data files and examples
- `configure.ac` - Autoconf configuration
- `Makefile.am` - Automake configuration

## Making Changes

1. **Edit org-mode files** in `/org` directory
2. **Write code in code blocks** with appropriate `#+begin_src` tags
3. **Document your code** thoroughly in the org-mode narrative
4. **Run make** to tangle and compile
5. **Run tests** with `make check`
6. **Never commit generated files** to the repository (only org files)

## Important Notes

- This is a library for scientific computing, precision and correctness are critical
- All changes should be accompanied by tests
- Documentation is as important as code (they're the same in org-mode!)
- Performance matters: consider both standard and HPC builds
- Thread safety: many functions are meant to be thread-safe
- The library is designed to never abort and never do I/O directly

## Resources

- Documentation: https://trex-coe.github.io/qmckl/index.html
- Source Repository: https://github.com/TREX-CoE/qmckl
- TREX Project: https://trex-coe.eu
