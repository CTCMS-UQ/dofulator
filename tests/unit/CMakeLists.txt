# Each unit test is a .c file.
# unit_test.c contains the test framework and compiles to a shared library

file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
file(GLOB UNIT_LIB_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/unit_test.c)
list(REMOVE_ITEM UNIT_TEST_SRCS ${UNIT_LIB_SRCS})

# libunittest.so
add_library(unittest SHARED ${UNIT_LIB_SRCS})
target_link_libraries(unittest PRIVATE dofulator m)
target_include_directories(unittest PUBLIC ${INCLUDE_DIRS})
target_compile_options(unittest PUBLIC "-Wall;-Wextra;-Wpedantic")

# unit test binaries
foreach(TEST_NAME IN ITEMS ${UNIT_TEST_SRCS})
  get_filename_component(TEST_BASENAME ${TEST_NAME} NAME_WLE)
  add_executable(${TEST_BASENAME} ${TEST_NAME})
  target_link_libraries(${TEST_BASENAME} PRIVATE unittest m)
  add_test(NAME unit::${TEST_BASENAME} COMMAND ${TEST_BASENAME})
endforeach()
