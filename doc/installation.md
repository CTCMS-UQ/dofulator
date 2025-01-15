# Installation

## For C projects and CLI

The library (`dofulator`) and CLI (`dof`) can be built directly with CMake:
```bash
mkdir build && cd build
cmake ..
cmake --build .
cmake --install . --component dofulator   # optional - install the library
cmake --install . --component dof         # optional - install the CLI
```
This will require CMake to be able to find a BLAS and LAPACK installation
(the MKL and OpenBLAS implementations have both been confirmed to work).
BLAS and LAPACK should be found automatically if they are installed,
but some environments may require additional configuration or specific flags to be set.
If you encounter errors with `FindBLAS`, this is likely the reason.

Alternatively, on all systems, the OpenBLAS implementation provided by the
[scipy-openblas32](https://pypi.org/project/scipy-openblas32/) Python package
is supported, and can be compiled against by setting `-DDOF_USE_SCIPY_OPENBLAS32=ON`.

#### Linux Users

When using system BLAS/LAPACK libraries, make sure you have installed the -dev/-devel
versions if your package manager uses split packages, and ensure that you have the
cblas and lapacke components installed.

In a HPC environment, this should simply require loading the desired BLAS module (both at build and run time).

#### Mac Users

On Mac OS, the Accelerate library is currently NOT supported.
Instead, if you'd rather not use scipy-openblas32, other BLAS implementations
(e.g. regular OpenBLAS) can be installed via Homebrew.
Since Mac doesn't add those packages to the system path by default, the
following flags (or similar) will likely be needed to properly locate the
installation:
```
cmake .. -DBLA_VENDOR=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/homebrew/Cellar/openblas/<VERSION_NUMBER>/ -DCMAKE_C_FLAGS='-I/opt/homebrew/Cellar/openblas/<VERSION_NUMBER/include'
```
Make sure to set the correct version number for the version you have installed.

#### Windows Users

The default MSVC C compiler does not fully comply with modern C standards, so it is not supported.
Instead, a simple alternative is to use the [Zig compile C](https://zighelp.org/chapter-4/#zig-cc-zig-c)
functionality by downloading the Zig compiler and setting `-DZIGLANG_PATH=C:\\path\to\ziglang\`
when first configuring CMake, or simply `-DCMAKE_C_COMPILER="zig;cc"` if the
`zig` binary is already findable by CMake (this is also an option on other platforms).

The other alternative is to use the gcc or clang compiler via [MinGW](https://www.mingw-w64.org/),
[Cygwin](https://cygwin.com/), or [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

### CMake Configuration

| Flag                      | Default   | Description                                           |
|---------------------------|-----------|-------------------------------------------------------|
| `DOF_LIBRARY_ONLY`        | `OFF`     | Build only the library (and not the CLI)      |
| `DOF_SHARED_LIB`          | `OFF`     | Build as a shared library (default is to build a static library) |
| `DOF_ATOM_TAG_TYPE`       | `int64_t` | Data type to store atom identifiers. Must be an integer type. |
| `DOF_ENABLE_PYTHON`       | `ON`      | Include Python targets and tests. |
| `DOF_USE_SCIPY_OPENBLAS32`| `OFF`     | Build against the scipy-openblas32 library, assumed to be installed in the current Python environment. `BLA_VENDOR` is ignored if `DOF_USE_SCIPY_OPENBLAS32` is set. |
| `BLA_VENDOR`              |           | Override the default BLAS vendor (see [FindBLAS](https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors)). |
| `FORCE_CBLAS_INCLUDE`     |           | Header file to include for CBLAS. The build script tries to auto-detect this, but if it fails for an obscure BLAS vendor it can be set manually here.|
| `FORCE_LAPACKE_INCLUDE`   |           | Header file to include for LAPACKE. The build script tries to auto-detect this, but if it fails for an obscure BLAS vendor it can be set manually here.|
| `ZIGLANG_PATH`            |           | Path to the Zig directory. Setting this will make CMake prioritise using the zig cc compiler. |
| `CMAKE_INSTALL_PREFIX`    |           | Directory to install library/CLI to with `cmake --install` |
| `ENABLE_LTO`              | `ON`      | Use link-time optimisations if supported |
| `ENABLE_TESTING`          | `OFF`     | Build tests (run with `ctest`) |
| `ENABLE_ASAN`             | `OFF`     | Build with address sanitization (for debugging purposes) |
| `ENABLE_UBSAN`            | `OFF`     | Build with undefined behaviour sanitization (for debugging purposes) |

### Testing

To run automated tests, from the `build` directory:
```bash
cmake . -DENABLE_TESTING=ON
cmake --build . --target tests
ctest
```

## For Python projects

With the repository downloaded, the package can be installed from the project root directory with
```bash
pip install .
```
or
```bash
pip install .[mdanalysis]
```
Alternatively, it can be installed directly from the git repository with
```bash
pip install "dofulator[mdanalysis]@git+https://github.com/CTCMS-UQ/dofulator.git"
```

As the C backend requires BLAS and LAPACK, these must be present on your system.
By default, the [scipy-openblas32](https://pypi.org/project/scipy-openblas32/)
package is installed and built against, but you can override this behaviour to
use your system-installed BLAS and LAPACK libraries with the flag
`--config-settings=cmake.define.DOF_USE_SCIPY_OPENBLAS32=OFF` (e.g. this may be
desirable for performance in a HPC environment).

MKL and OpenBLAS have both been tested to work. Other implementations should
also work, but may require some additional configuration (please submit an
issue if you encounter problems).
If the default implementation found is not the one you want, a particular
system BLAS implementation can be chosen on installation, for example to
choose MKL:
```bash
pip install . --config-settings=cmake.define.BLA_VENDOR=MKL --config-settings=cmake.define.DOF_USE_SCIPY_OPENBLAS32=OFF
```
See the [FindBLAS](https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors)
documentation for details of other vendor options.

If your BLAS implementation uses odd names for the `cblas.h` or `lapacke.h` headers, you
can set these directly with the flags `--config-settings=cmake.define.FORCE_CBLAS_INCLUDE=my_cblas.h`
and `--config-settings=cmake.define.FORCE_LAPECKE_INCLUDE=my_lapacke.h`
(replacing `my_cblas.h` and `my_lapacke.h` with the appropriate header names).

On Windows, to avoid incompatibilities with the MSVC compiler, the Zig compiler bundled in
the [ziglang Python package](https://pypi.org/project/ziglang/) is used by default.
Other compilers can be selected by setting `--config-settings=cmake.define.ZIGLANG_PATH=""`
along with the other appropriate flags for your compiler of choice.
