# Installation

## For C projects and CLI

The library (`dofulator`) and CLI (`dof`) can be built directly with CMake:
```bash
$ mkdir build && cd build
$ cmake ..
$ cmake --build .
$ cmake --install . --component dofulator   # optional - install the library
$ cmake --install . --component dof         # optional - install the CLI
```
This will require CMake to be able to find a BLAS and LAPACK installation
(the MKL and OpenBLAS implementations have both been confirmed to work).
BLAS and LAPACK should be found automatically if they are installed,
but some environments may require additional configuration or specific flags to be set.
If you encounter errors with `FindBLAS`, this is likely the reason.

### CMake Configuration

| Flag                    | Default   | Description                                           |
|-------------------------|-----------|-------------------------------------------------------|
| `DOF_LIBRARY_ONLY`      | `OFF`     | Build only the library (and not the CLI)      |
| `DOF_SHARED_LIB`        | `OFF`     | Build as a shared library (default is to build a static library) |
| `DOF_ATOM_TAG_TYPE`     | `int64_t` | Data type to store atom identifiers. Must be an integer type. |
| `DOF_ENABLE_PYTHON`     | `ON`      | Include Python targets and tests. |
| `BLA_VENDOR`            |           | Override the default BLAS vendor (see [FindBLAS](https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors)). |
| `FORCE_CBLAS_INCLUDE`   |           | Header file to include for CBLAS. The build script tries to auto-detect this, but if it fails for an obscure BLAS vendor it can be set manually here.|
| `FORCE_LAPACKE_INCLUDE` |           | Header file to include for LAPACKE. The build script tries to auto-detect this, but if it fails for an obscure BLAS vendor it can be set manually here.|
| `ENABLE_LTO`            | `ON`      | Use link-time optimisations if supported |
| `CMAKE_INSTALL_PREFIX`  |           | Directory to install library/CLI to with `cmake --install` |
| `ENABLE_TESTING`        | `OFF`     | Build tests (run with `ctest`) |
| `ENABLE_ASAN`           | `OFF`     | Build with address sanitization (for debugging purposes) |
| `ENABLE_UBSAN`          | `OFF`     | Build with undefined behaviour sanitization (for debugging purposes) |

### Testing

To run automated tests, from the `build` directory:
```bash
$ cmake . -DENABLE_TESTING=ON
$ cmake --build . --target tests
$ ctest
```

## For Python projects

With the repository downloaded, the package can be installed from the project root directory with
```bash
pip install .
```
Alternatively, it can be installed directly from the git repository with
```bash
pip install git+https://github.com/CTCMS-UQ/dofulator
```

As the C backend requires BLAS and LAPACK, these must be present on your system.
You can install them with your system package manager, or via other package
managers like conda.
MKL and OpenBLAS have both been tested to work. Other implementations should
also work, but may require some additional configuration (please submit an
issue if you encounter problems).
A particular BLAS implementation can be chosen on installation, for example to
choose OpenBLAS:
```bash
pip install . --config-settings=cmake.define.BLA_VENDOR=OpenBLAS
```
See the [FindBLAS](https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors)
documentation for details of other vendor options.

If your BLAS implementation uses odd names for the `cblas.h` or `lapacke.h` headers, you
can set these directly with the flags `--config-settings=cmake.define.FORCE_CBLAS_INCLUDE=my_cblas.h`
and `--config-settings=cmake.define.FORCE_LAPECKE_INCLUDE=my_lapacke.h`
(replacing `my_cblas.h` and `my_lapacke.h` with the appropriate header names).
