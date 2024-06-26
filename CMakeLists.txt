cmake_minimum_required(VERSION 3.13)
project(dofulator C)

if (NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
endif()
if (CMAKE_C_STANDARD LESS 11)
  message(FATAL_ERROR "C standard must be set to at least 11")
endif()

set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF CACHE BOOL "Use compiler extensions")

set(ENABLE_LTO ON CACHE BOOL "Enable link-time optimisation if supported")

set(ENABLE_TESTING OFF CACHE BOOL "Compile unit and integration tests")
set(ENABLE_ASAN OFF CACHE BOOL "Add -fsanitize=address -fno-omit-frame-pointer")
set(ENABLE_UBSAN OFF CACHE BOOL "Add -fsanitize=undefined")

set(DOF_LIBRARY_ONLY OFF CACHE BOOL "Only build the library")
set(DOF_SHARED_LIB OFF CACHE BOOL "Build library as shared. Default is static.")
set(DOF_ATOM_TAG_TYPE "int64_t" CACHE STRING "Integer data type for atom identifiers. Default int64_t required for compatibility with python wrapper.")

set(FORCE_CBLAS_INCLUDE "" CACHE STRING "Override name of cblas header to include")
set(FORCE_LAPACKE_INCLUDE "" CACHE STRING "Override name of lapacke header to include")

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
target_include_directories(dofulator PRIVATE ${INCLUDE_DIRS})

# Warnings
include(CheckCCompilerFlag)
foreach(FLAG IN ITEMS -Wall -Wextra -Wpedantic)
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
set(BLA_STATIC ON)
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
  set_cblas_lapacke_headers("openblas/cblas.h" "openblas/lapacke.h")
  target_compile_definitions(dofulator PRIVATE BLA_VENDOR_OPENBLAS)
elseif (DETECTED_VENDOR STREQUAL MKL)
  find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
  set_cblas_lapacke_headers("mkl_cblas.h" "mkl_lapacke.h")
  target_compile_options(dofulator PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_compile_definitions(dofulator PRIVATE MKL_EXTENSIONS)
  target_include_directories(dofulator PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(dofulator PUBLIC $<LINK_ONLY:MKL::MKL>)
else()
  message(WARNING
    "Unsupported BLAS/LAPACK vendor - may not be able to find headers."
    "If compilation fails, try pointing INCLUDE_DIRECTORIES to the include directory of the BLAS/LAPACK implementations"
    "and/or setting FORCE_CBLAS_INCLUDE and FORCE_LAPACKE_INCLUDE to the correct header file names for the implementation."
    "Alternatively, you can set BLA_VENDOR to use a different BLAS/LAPACK implementation."
  )
  set_cblas_lapacke_headers("cblas.h" "lapacke.h")
endif()

configure_file(${SRC_DIR}/cblas_lapacke.in.h ${CMAKE_BINARY_DIR}/cblas_lapacke.h @ONLY USE_SOURCE_PERMISSIONS)
target_include_directories(dofulator PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(dofulator PRIVATE m ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_options(dofulator PRIVATE ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})


################################################################################
# Address/UB sanitization
if (ENABLE_ASAN)
  target_link_libraries(dofulator PUBLIC -fsanitize=address -fno-omit-frame-pointer)
endif()

if (ENABLE_UBSAN)
  target_link_libraries(dofulator PUBLIC -fsanitize=undefined)
endif()


################################################################################
# dof binary
if (NOT DOF_LIBRARY_ONLY)
  file(GLOB EXE_SRCS ${SRC_DIR}/frontend/*.c)
  add_executable(dof ${EXE_SRCS})
  target_include_directories(dof PRIVATE ${SRC_DIR})
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
# Set up testing
if (ENABLE_TESTING)
  message(STATUS "Testing: enabled")
  if (ENABLE_ASAN)
    message(STATUS "Address sanitizer in unit tests: enabled")
  else()
    message(STATUS "Address sanitizer in unit tests: disabled")
  endif()
  include(CTest)
  add_subdirectory(tests)
else()
  message(STATUS "Testing: disabled")
endif()


################################################################################
# Optional python package
find_package(Python COMPONENTS Interpreter Development.Module NumPy)
if (Python_FOUND)
  message(STATUS "Python package target available")
  add_subdirectory(python)
endif()
