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
set(ENABLE_TESTING ON CACHE BOOL "Compile unit and integration tests")
set(TEST_ASAN OFF CACHE BOOL "Add -fsanitize=address to unit and integration test binaries")

set(INCLUDE_DIRS ${INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src)
file(GLOB_RECURSE LIB_SRCS ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB EXE_SRCS ${PROJECT_SOURCE_DIR}/src/main.c)
list(REMOVE_ITEM LIB_SRCS ${EXE_SRCS})


add_library(dofulator STATIC ${LIB_SRCS})
target_include_directories(dofulator INTERFACE ${INCLUDE_DIRS})
target_compile_options(dofulator PRIVATE "-Wall;-Wextra;-Wpedantic")

install(TARGETS dofulator EXPORT dofulator_Targets DESTINATION lib)
install(FILES ${PROJECT_SOURCE_DIR}/src/dofulator.h DESTINATION include)

find_package(LAPACK)
find_package(BLAS)

target_link_libraries(dofulator PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


add_executable(dof ${EXE_SRCS})
target_link_libraries(dof PRIVATE dofulator m)

if (ENABLE_TESTING)
  include(CTest)

  file(GLOB UNIT_TEST_SRCS ${PROJECT_SOURCE_DIR}/tests/unit/*.c)
  file(GLOB UNIT_LIB_SRCS ${PROJECT_SOURCE_DIR}/tests/unit/unit_test.c)
  list(REMOVE_ITEM UNIT_TEST_SRCS ${UNIT_LIB_SRCS})

  add_library(unittest SHARED ${UNIT_LIB_SRCS})
  target_link_libraries(unittest PRIVATE dofulator m)
  target_include_directories(unittest PUBLIC ${INCLUDE_DIRS})
  target_compile_options(unittest PUBLIC "-Wall;-Wextra;-Wpedantic")
  if (TEST_ASAN)
    target_link_libraries(unittest PUBLIC -fsanitize=address)
  endif()

  foreach(TEST_NAME IN ITEMS ${UNIT_TEST_SRCS})
    get_filename_component(TEST_BASENAME ${TEST_NAME} NAME_WLE)
    add_executable(${TEST_BASENAME} ${TEST_NAME})
    target_link_libraries(${TEST_BASENAME} PRIVATE unittest)
    add_test(NAME ${TEST_BASENAME} COMMAND ${TEST_BASENAME})
  endforeach()

endif()
