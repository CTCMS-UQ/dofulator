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
set(ENABLE_TESTING ON CACHE BOOL "Compile unit and integration tests")
set(TEST_ASAN OFF CACHE BOOL "Add -fsanitize=address to unit and integration test binaries")

set(INCLUDE_DIRS ${INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src)
file(GLOB_RECURSE LIB_SRCS ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB EXE_SRCS ${PROJECT_SOURCE_DIR}/src/main.c)
list(REMOVE_ITEM LIB_SRCS ${EXE_SRCS})


# libdofulator.a
add_library(dofulator STATIC ${LIB_SRCS})
target_include_directories(dofulator INTERFACE ${INCLUDE_DIRS})
target_compile_options(dofulator PRIVATE "-Wall;-Wextra;-Wpedantic")

install(TARGETS dofulator EXPORT dofulator_Targets DESTINATION lib)
install(FILES ${PROJECT_SOURCE_DIR}/src/dofulator.h DESTINATION include)

find_package(LAPACK)
find_package(BLAS)

target_link_libraries(dofulator PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


# dof binary
add_executable(dof ${EXE_SRCS})
target_link_libraries(dof PRIVATE dofulator m)


# Turn on link time optimisation if supported
if (ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)

  if(LTO_SUPPORTED)
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET dofulator PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET dof PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO / LTO not supported: <${LTO_ERROR}>")
  endif()
else()
  message(STATUS "IPO / LTO disabled")
endif()

# Set up testing
if (ENABLE_TESTING)
  include(CTest)
  add_subdirectory(tests)
endif()
