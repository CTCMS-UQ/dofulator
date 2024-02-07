cmake_minimum_required(VERSION 3.10)
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
set(TEST_ASAN OFF CACHE BOOL "Add -fsanitize=address to unit test binaries")
set(DOF_LIBRARY_ONLY OFF CACHE BOOL "Only build the library")

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${SRC_DIR})
file(GLOB LIB_SRCS ${SRC_DIR}/*.c)
file(GLOB EXE_SRCS ${SRC_DIR}/frontend/*.c)


# libdofulator.a
add_library(dofulator STATIC ${LIB_SRCS})
target_include_directories(dofulator PRIVATE ${INCLUDE_DIRS})
target_compile_options(dofulator PRIVATE "-Wall;-Wextra;-Wpedantic")

install(TARGETS dofulator EXPORT dofulator_Targets DESTINATION lib)
install(FILES ${SRC_DIR}/dofulator.h ${SRC_DIR}/atom_list.h ${SRC_DIR}/bond_list.h DESTINATION include)

set(BLA_STATIC ON)
find_package(LAPACK)
find_package(BLAS)

target_link_libraries(dofulator PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_options(dofulator PRIVATE ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})


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
  if (TEST_ASAN)
    message(STATUS "Address sanitizer in unit tests: enabled")
  else()
    message(STATUS "Address sanitizer in unit tests: disabled")
  endif()
  include(CTest)
  add_subdirectory(tests)
else()
  message(STATUS "Testing: disabled")
endif()
