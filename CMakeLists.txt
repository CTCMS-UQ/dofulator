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

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIRS ${INCLUDE_DIRECTORIES} ${SRC_DIR})
# file(GLOB LIB_SRCS ${SRC_DIR}/*.c)
file(GLOB LIB_SRCS ${SRC_DIR}/dofulator.c)
file(GLOB EXE_SRCS ${SRC_DIR}/frontend/*.c)


# libdofulator.a
add_library(dofulator STATIC ${LIB_SRCS})
target_include_directories(dofulator PRIVATE ${INCLUDE_DIRS})
target_compile_options(dofulator PRIVATE "-Wall;-Wextra;-Wpedantic")

install(TARGETS dofulator EXPORT dofulator_Targets DESTINATION lib)
# install(FILES ${SRC_DIR}/dofulator.h ${SRC_DIR}/atom_list.h ${SRC_DIR}/bond_list.h DESTINATION include)
install(FILES ${SRC_DIR}/dofulator.h DESTINATION include)

# Set up cblas and lapacke includes
# Default vendor can be overridden by setting BLA_VENDOR
# OpenBLAS seems to beat MKL in speed for all test cases
set(BLA_STATIC ON)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

if (BLA_VENDOR STREQUAL OpenBLAS)
  set(INCLUDE_CBLAS "#include <openblas/cblas.h>")
  set(INCLUDE_LAPACKE "#include <openblas/lapacke.h>")
elseif (BLA_VENDOR MATCHES "[(Intel)(Intel10_32)(Intel10_64lp)(Intel10_64lp_seq)(Intel10_64ilp)(Intel10_64ilp_seq)(Intel10_64_dyn)]")
  find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
  set(INCLUDE_CBLAS "#include <mkl_cblas.h>")
  set(INCLUDE_LAPACKE "#include <mkl_lapacke.h>")
  target_compile_options(dofulator PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_compile_definitions(dofulator PRIVATE MKL_EXTENSIONS)
  target_include_directories(dofulator PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(dofulator PUBLIC $<LINK_ONLY:MKL::MKL>)
else()
  message(WARNING "Unsupported blas vendor - may not be able to find the header."
    "If compilation fails, try pointing INCLUDE_DIRECTORIES to the include directory of the blas implementation.")
  set(INCLUDE_CBLAS "#include <cblas.h>")
  set(INCLUDE_LAPACKE "#include <lapacke.h>")
endif()

configure_file(${SRC_DIR}/cblas_lapacke.h.in ${CMAKE_BINARY_DIR}/cblas_lapacke.h @ONLY USE_SOURCE_PERMISSIONS)
target_include_directories(dofulator PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(dofulator PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_options(dofulator PRIVATE ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
if (ENABLE_ASAN)
  target_link_libraries(dofulator PUBLIC -fsanitize=address -fno-omit-frame-pointer)
endif()
if (ENABLE_UBSAN)
  target_link_libraries(dofulator PUBLIC -fsanitize=undefined)
endif()


# dof binary
if (NOT DOF_LIBRARY_ONLY)
  add_executable(dof ${EXE_SRCS})
  target_include_directories(dof PRIVATE ${SRC_DIR})
  target_link_libraries(dof PRIVATE dofulator m)
endif()


# Turn on link time optimisation if supported
if (ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)

  if(LTO_SUPPORTED)
    message(STATUS "IPO / LTO: enabled")
    set_property(TARGET dofulator PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET dof PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO / LTO: not supported - <${LTO_ERROR}>")
  endif()
else()
  message(STATUS "IPO / LTO: disabled")
endif()

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
