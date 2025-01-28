cmake_minimum_required(VERSION 3.25)

################################################################################

set(ENABLE_LTO ON CACHE BOOL "Enable link-time optimisation if supported")

set(ENABLE_TESTING OFF CACHE BOOL "Compile unit and integration tests")
set(ENABLE_ASAN OFF CACHE BOOL "Add -fsanitize=address -fno-omit-frame-pointer")
set(ENABLE_UBSAN OFF CACHE BOOL "Add -fsanitize=undefined")

set(DOF_LIBRARY_ONLY OFF CACHE BOOL "Only build the library")
set(DOF_SHARED_LIB OFF CACHE BOOL "Build library as shared. Default is static.")
set(DOF_ATOM_TAG_TYPE "int64_t" CACHE STRING "Integer data type for atom identifiers. Default int64_t required for compatibility with python wrapper.")
set(DOF_ENABLE_PYTHON ON CACHE BOOL "Enable building/testing of the python interface")

set(DOF_USE_SCIPY_OPENBLAS32 OFF CACHE BOOL "Build against scipy-openblas32 BLAS/LAPACK library")
set(FORCE_CBLAS_INCLUDE "" CACHE STRING "Override name of cblas header to include")
set(FORCE_LAPACKE_INCLUDE "" CACHE STRING "Override name of lapacke header to include")

################################################################################

# The msvc C compiler does not support foo[restrict N] syntax or pointer to VLA types.
# For portability, the zig cc compiler can be used either from wherever it's installed
# in the system or from the ziglang python package.
# Setting the ZIGLANG_PATH variable achieves the former, or
# DOF_USE_PY_ZIGLANG=ON achieves the latter.
#
# Need to do this before calling project()

set(ZIGLANG_PATH "" CACHE STRING "Path to the Zig directory. Setting this will make CMake prioritise using the zig cc compiler.")
set(DOF_USE_PY_ZIGLANG OFF CACHE BOOL "Build with the zig C compiler bundled with the ziglang Python package")
set(PATCH_ZIGLANG_0_13_0_STACK_PROBE OFF CACHE BOOL "Automatically patch ziglang to avoid a linkage bug on Windows.")
set(TRY_USING_ZIG OFF CACHE BOOL "Prioritise compiling with zig cc")

# Look for zig in ZIGLANG_PATH if defined
if (NOT ${ZIGLANG_PATH} STREQUAL "")
  list(PREPEND CMAKE_PREFIX_PATH ${ZIGLANG_PATH})
  set(TRY_USING_ZIG ON CACHE BOOL "" FORCE)
endif()

# Look for ziglang python package
if (${DOF_USE_PY_ZIGLANG})
  set(TRY_USING_ZIG ON CACHE BOOL "" FORCE)
  find_package(Python COMPONENTS Interpreter)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import ziglang" RESULT_VARIABLE CHECK_ZIGLANG OUTPUT_QUIET ERROR_QUIET)
  if (CHECK_ZIGLANG EQUAL 0)
    message(STATUS "Found ziglang python package. Compiling with zig cc")
    list(PREPEND CMAKE_PREFIX_PATH ${Python_SITELIB}/ziglang)
  else()
    message(WARNING "Could not find ziglang python package! Falling back to default compiler.")
  endif()
endif()

if (TRY_USING_ZIG)
  # Under scikit-build, site-packages lib for the isolated build environment should
  # be in the cmake path, so ziglang/zig should match.
  # Fallback to searching for zig in case not being run via scikit-build
  find_program(ZIG_PATH ziglang/zig)
  if (ZIG_PATH STREQUAL "ZIG_PATH-NOTFOUND")
    find_program(ZIG_PATH zig)
  endif()
  if (ZIG_PATH STREQUAL "ZIG_PATH-NOTFOUND")
    message(WARNING "Could not find the zig compiler! Falling back to default.")
  else()
    # Wrappers for ar and ranlib since cmake doesn't support definition of extra arguments
    macro(ZIGCMD CMD CMAKE_VAR)
      if (WIN32)
        write_file(zig${CMD}.cmd "${ZIG_PATH} ${CMD} %*")
        set(${CMAKE_VAR} ${CMAKE_BINARY_DIR}/zig${CMD}.cmd)
      else()
        write_file(zig${CMD} "ZIG_LOCAL_CACHE_DIR=${CMAKE_BINARY_DIR}/zig_cache ZIG_GLOBAL_CACHE_DIR=${CMAKE_BINARY_DIR}/zig_cache ${ZIG_PATH} ${CMD} \"$@\"")
        file(CHMOD zig${CMD} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)
        set(${CMAKE_VAR} ${CMAKE_BINARY_DIR}/zig${CMD})
      endif()
    endmacro()
    ZIGCMD(cc CMAKE_C_COMPILER)
    ZIGCMD(ar CMAKE_C_COMPILER_AR)
    ZIGCMD(ranlib CMAKE_C_COMPILER_RANLIB)

    # The latest ziglang package (0.13.0.post1) has a bug for linkage on windows which will
    # be fixed in the next release. Apply the patch manually for now if enabled.
    if (WIN32)
      execute_process(COMMAND ${ZIG_PATH} version OUTPUT_VARIABLE ZIG_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
      if (ZIG_VERSION VERSION_EQUAL 0.13.0)
        if (PATCH_ZIGLANG_0_13_0_STACK_PROBE)
          # Get path to zig and patch stack_probe.zig
          cmake_path(REMOVE_FILENAME ZIG_PATH OUTPUT_VARIABLE ZIGLANG_DIR)
          message(WARNING "Patching zig's ${ZIGLANG_DIR}/lib/compiler_rt/stack_probe.zig file to https://github.com/rsepassi/zig/blob/ca589fd41e37fd1dbd863eea33da9365ac7a3ac5/lib/compiler_rt/stack_probe.zig to avoid a bug. See https://github.com/ziglang/zig/pull/20138 and https://github.com/ziglang/zig/issues/15107 for details. The old version of the file will be stored as ${CMAKE_BINARY_DIR}\\stack_probe.zig.old")
          file(DOWNLOAD https://raw.githubusercontent.com/rsepassi/zig/ca589fd41e37fd1dbd863eea33da9365ac7a3ac5/lib/compiler_rt/stack_probe.zig stack_probe.zig EXPECTED_MD5;aa82aa00300cc51aa0aa9f702c43e903 SHOW_PROGRESS)
          execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${ZIGLANG_DIR}/lib/compiler_rt/stack_probe.zig stack_probe.zig.old)
          execute_process(COMMAND ${CMAKE_COMMAND} -E copy stack_probe.zig ${ZIGLANG_DIR}/lib/compiler_rt/stack_probe.zig)
        else()
          message(WARNING "Zig version ${ZIG_VERSION} detected which has a linkage bug on Windows (https://github.com/ziglang/zig/pull/20138). If compilation fails, please use a later version, patch as per the github link, or enable automatic patching by setting PATCH_ZIGLANG_0_13_0_STACK_PROBE=ON.")
        endif()
      elseif(ZIG_VERSION VERSION_LESS 0.13.0)
        message(WARNING "Zig version ${ZIG_VERSION} detected which has a linkage bug on Windows (https://github.com/ziglang/zig/pull/20138). If compilation fails, please use a later version or patch as per the github link, or use version 0.13.0 with PATCH_ZIGLANG_0_13_0_STACK_PROBE=ON to enable automated patching.")
      endif()
    endif()
  endif()
endif()

################################################################################

project(dofulator C)

if (NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
endif()
if (CMAKE_C_STANDARD LESS 11)
  message(FATAL_ERROR "C standard must be set to at least 11")
endif()

set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF CACHE BOOL "Use compiler extensions")

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIRS ${INCLUDE_DIRECTORIES} ${SRC_DIR})

################################################################################
# dofulator core library
file(GLOB LIB_SRCS ${SRC_DIR}/dofulator.c)
if (DOF_SHARED_LIB)
  message(STATUS "Building dofulator backend as shared library")
  add_library(dofulator SHARED ${LIB_SRCS})
else()
  message(STATUS "Building dofulator backend as static library")
  add_library(dofulator STATIC ${LIB_SRCS})
endif()
target_include_directories(dofulator BEFORE PRIVATE ${INCLUDE_DIRS})

# Warnings
include(CheckCCompilerFlag)
foreach(FLAG IN ITEMS -Wall -Wextra -Wpedantic -fPIC)
  check_c_compiler_flag(${FLAG} HAVE_FLAG)
  if (HAVE_FLAG)
    target_compile_options(dofulator PRIVATE ${FLAG})
  endif()
endforeach()

if (NOT ${DOF_ATOM_TAG_TYPE} EQUAL "")
  message(STATUS "AtomTag type: ${DOF_ATOM_TAG_TYPE}")
  target_compile_definitions(dofulator PRIVATE AtomTag=${DOF_ATOM_TAG_TYPE})
else()
  message(STATUS "AtomTag type: int64_t")
endif()

install(TARGETS dofulator EXPORT dofulator_Targets DESTINATION lib)
install(FILES ${SRC_DIR}/dofulator.h DESTINATION include)

################################################################################
# Set up cblas and lapacke includes/linkage
# Default vendor can be overridden by setting BLA_VENDOR
# OpenBLAS seems to marginally beat MKL in speed for all test cases

# Force use of scipy-openblas if specified
if (${DOF_USE_SCIPY_OPENBLAS32})
  find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
  set(INCLUDE_CBLAS "cblas.h")
  set(INCLUDE_LAPACKE "lapacke.h")
  # Get details from python package
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import scipy_openblas32 as spob;print(f'{spob.get_library()};{spob.get_lib_dir()};{spob.get_include_dir()}',end='')" OUTPUT_VARIABLE SCIPY_OPENBLAS_DETAILS)
  list(POP_BACK SCIPY_OPENBLAS_DETAILS SCIPY_OPENBLAS_INCLUDE)
  list(POP_BACK SCIPY_OPENBLAS_DETAILS SCIPY_OPENBLAS_LIB)
  list(POP_BACK SCIPY_OPENBLAS_DETAILS BLAS_LIBRARIES)
  target_include_directories(dofulator PUBLIC ${SCIPY_OPENBLAS_INCLUDE})
  target_link_directories(dofulator PUBLIC ${SCIPY_OPENBLAS_LIB})
  target_compile_definitions(dofulator PUBLIC SCIPY_OPENBLAS)
  message(STATUS "Using scipy-openblas32: ${SCIPY_OPENBLAS_LIB}")
else()

  find_package(LAPACK REQUIRED)
  find_package(BLAS REQUIRED)

# Detect found vendor
  set(DETECTED_VENDOR "Unknown")
  foreach (BLIB ${BLAS_LIBRARIES})
    if (${BLIB} MATCHES openblas)
      set(DETECTED_VENDOR "OpenBLAS")
      break()
    elseif (${BLIB} MATCHES mkl)
      set(DETECTED_VENDOR "MKL")
      break()
    endif()
  endforeach()
  if (DETECTED_VENDOR STREQUAL "Unknown" AND EXISTS ${BLAS_Accelerate_LIBRARY})
    set(DETECTED_VENDOR "Accelerate")
  endif()
  message(STATUS "Detected BLAS vendor: ${DETECTED_VENDOR}")

  macro(set_cblas_lapacke_headers CBLAS LAPACKE)
    if (FORCE_CBLAS_INCLUDE STREQUAL "")
      set(INCLUDE_CBLAS ${CBLAS})
    else()
      set(INCLUDE_CBLAS ${FORCE_CBLAS_INCLUDE})
    endif()
    if (FORCE_LAPACKE_INCLUDE STREQUAL "")
      set(INCLUDE_LAPACKE ${LAPACKE})
    else()
      set(INCLUDE_LAPACKE ${FORCE_LAPACKE_INCLUDE})
    endif()
  endmacro()

  if (DETECTED_VENDOR STREQUAL OpenBLAS)
    # Some package managers put openblas headers under a prefix, so detect that
    set_cblas_lapacke_headers("openblas/cblas.h" "openblas/lapacke.h")
    try_compile(OPENBLAS_PREFIX_INCLUDE
      SOURCE_FROM_CONTENT
        openblas_test.c
        "#include <${INCLUDE_CBLAS}>\n#include <${INCLUDE_LAPACKE}>\nint main(void) {}"
      NO_CACHE
      LINK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
      LINK_OPTIONS ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}
    )
    if (NOT OPENBLAS_PREFIX_INCLUDE)
      set_cblas_lapacke_headers("cblas.h" "lapacke.h")
    endif()
    target_compile_definitions(dofulator PRIVATE BLA_VENDOR_OPENBLAS)
  elseif (DETECTED_VENDOR STREQUAL MKL)
    find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
    set_cblas_lapacke_headers("mkl_cblas.h" "mkl_lapacke.h")
    target_compile_options(dofulator PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_compile_definitions(dofulator PRIVATE MKL_EXTENSIONS)
    target_include_directories(dofulator PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(dofulator PUBLIC $<LINK_ONLY:MKL::MKL>)
  elseif(DETECTED_VENDOR STREQUAL Accelerate)
    set_cblas_lapacke_headers("Accelerate/Accelerate.h" "Accelerate/Accelerate.h")
    target_compile_definitions(dofulator PRIVATE -DACCELERATE_NEW_LAPACK -Dlapack_int=__LAPACK_int)
  else()
    message(WARNING
      "Unsupported BLAS/LAPACK vendor - may not be able to find headers."
      "If compilation fails, try pointing INCLUDE_DIRECTORIES to the include directory of the BLAS/LAPACK implementations"
      "and/or setting FORCE_CBLAS_INCLUDE and FORCE_LAPACKE_INCLUDE to the correct header file names for the implementation."
      "Alternatively, you can set BLA_VENDOR to use a different BLAS/LAPACK implementation."
    )
    set_cblas_lapacke_headers("cblas.h" "lapacke.h")
  endif()
endif()

configure_file(${SRC_DIR}/cblas_lapacke.in.h ${CMAKE_BINARY_DIR}/cblas_lapacke.h @ONLY USE_SOURCE_PERMISSIONS)
target_include_directories(dofulator BEFORE PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(dofulator PRIVATE m ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_options(dofulator PRIVATE ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})


################################################################################
# Check for [[nodiscard]] compiler extension
try_compile(SUPPORT_NODISCARD
  SOURCE_FROM_CONTENT
    nodiscard_test.c
    "__extension__ typedef enum [[nodiscard]] Test { T1, T2 } Test; int main(void) { return T1; }"
  NO_CACHE
)
if (SUPPORT_NODISCARD)
  target_compile_definitions(dofulator PRIVATE HAS_NODISCARD)
endif()


################################################################################
# Address/UB sanitization
if (ENABLE_ASAN)
  target_link_libraries(dofulator PUBLIC -fsanitize=address -fno-omit-frame-pointer)
endif()

if (ENABLE_UBSAN)
  target_link_libraries(dofulator PUBLIC -fsanitize=undefined)
endif()

# zig cc enables asan/ubsan by default, so turn off if disabled
get_filename_component(COMPILER_BASENAME ${CMAKE_C_COMPILER} NAME)
if (${COMPILER_BASENAME} STREQUAL "zig")
  set(NO_SANITIZE "")
  if (NOT ENABLE_UBSAN)
    list(APPEND NO_SANITIZE "undefined")
  endif()
  if (NOT ENABLE_ASAN)
    list(APPEND NO_SANITIZE "address")
  endif()
  list(JOIN NO_SANITIZE "," NOSAN)
  if (NOT ${NOSAN} STREQUAL "")
    list(APPEND CMAKE_C_FLAGS "-fno-sanitize=${NOSAN}")
  endif()
endif()

################################################################################
# dof binary
if (NOT DOF_LIBRARY_ONLY)
  file(GLOB EXE_SRCS ${SRC_DIR}/frontend/*.c)
  add_executable(dof ${EXE_SRCS})
  target_include_directories(dof BEFORE PRIVATE ${SRC_DIR})
  target_link_libraries(dof PRIVATE dofulator m)
endif()


################################################################################
# Turn on link time optimisation if supported
if (ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)

  if(LTO_SUPPORTED)
    message(STATUS "IPO / LTO: enabled")
    set_property(TARGET dofulator PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    if (NOT DOF_LIBRARY_ONLY)
        set_property(TARGET dof PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
  else()
    message(STATUS "IPO / LTO: not supported - <${LTO_ERROR}>")
  endif()
else()
  message(STATUS "IPO / LTO: disabled")
endif()


################################################################################
# Optional python package
if (DOF_ENABLE_PYTHON)
  find_package(Python COMPONENTS Interpreter Development.Module NumPy)
  if (Python_FOUND)
    message(STATUS "Python package target available")
    add_subdirectory(python)
  endif()
endif()

################################################################################
# Set up testing
if (ENABLE_TESTING)
  message(STATUS "Testing: enabled")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")
  if (ENABLE_ASAN)
    message(STATUS "Address sanitizer in unit tests: enabled")
  else()
    message(STATUS "Address sanitizer in unit tests: disabled")
  endif()
  include(CTest)
  add_subdirectory(tests)
  message(STATUS "NOTE: Make sure to build the 'tests' target before running ctest!")
  list(POP_BACK CMAKE_MESSAGE_INDENT)
else()
  message(STATUS "Testing: disabled")
endif()

