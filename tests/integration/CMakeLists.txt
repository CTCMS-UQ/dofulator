# Each test is a .c file.
# integration_test.c contains the test framework and compiles to a shared library

file(GLOB INT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB TEST_LIB_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/integration_test.c)
list(REMOVE_ITEM INT_TEST_SRCS ${TEST_LIB_SRCS})

# shared library
add_library(inttest SHARED ${TEST_LIB_SRCS})
target_link_libraries(inttest PRIVATE dofulator m)
target_include_directories(inttest PUBLIC ${INCLUDE_DIRS})
foreach(FLAG IN ITEMS -Wall -Wextra -Wpedantic)
  check_c_compiler_flag(${FLAG} HAVE_FLAG)
  if (HAVE_FLAG)
    target_compile_options(inttest PUBLIC ${FLAG})
  endif()
endforeach()

# integration test binaries
foreach(TEST_NAME IN ITEMS ${INT_TEST_SRCS})
  get_filename_component(TEST_BASENAME ${TEST_NAME} NAME_WLE)
  add_executable(${TEST_BASENAME} ${TEST_NAME})
  target_link_libraries(${TEST_BASENAME} PRIVATE inttest m)
  add_test(NAME integration::${TEST_BASENAME} COMMAND ${TEST_BASENAME})
  add_dependencies(test_deps ${TEST_BASENAME})
endforeach()
