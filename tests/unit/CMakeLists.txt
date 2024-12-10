# Each test is a .c file.

file(GLOB UNIT_TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

foreach(TEST_NAME IN ITEMS ${UNIT_TEST_SRCS})
  get_filename_component(TEST_BASENAME ${TEST_NAME} NAME_WLE)
  add_executable(${TEST_BASENAME} ${TEST_NAME})
  target_link_libraries(${TEST_BASENAME} PRIVATE dofulator m)
  target_include_directories(${TEST_BASENAME} PRIVATE ${INCLUDE_DIRS})

  foreach(FLAG IN ITEMS -Wall -Wextra -Wpedantic)
    check_c_compiler_flag(${FLAG} HAVE_FLAG)
    if (HAVE_FLAG)
      target_compile_options(${TEST_BASENAME} PRIVATE ${FLAG})
    endif()
  endforeach()

  add_test(NAME unit::${TEST_BASENAME} COMMAND ${TEST_BASENAME})
  add_dependencies(test_deps ${TEST_BASENAME})
endforeach()
