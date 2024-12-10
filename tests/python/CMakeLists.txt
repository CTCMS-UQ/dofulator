

# Check for usable pytest in python environment, or fall back to Pytest CMake module
execute_process(COMMAND ${Python_EXECUTABLE} -m pytest -V RESULT_VARIABLE PYTEST_SUCCESS OUTPUT_VARIABLE PYTEST_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT "${PYTEST_SUCCESS}" STREQUAL "0")
  find_package(Pytest)
  if (Pytest_FOUND)
    set(PYTEST_CMD Pytest::Pytest)
  endif()
else()
  set(PYTEST_CMD Python::Interpreter -m pytest)
  message(STATUS "Found installed module ${PYTEST_VERSION}")
endif()

if (PYTEST_CMD)
  add_test(
    NAME python_tests
    COMMAND ${PYTEST_CMD} ${CMAKE_CURRENT_SOURCE_DIR}
  )

  # Set up fixture to copy c_dofulator library into source tree/clean it up for pytest
  add_dependencies(test_deps c_dofulator)
  set(PYTHON_SRC_DIR ${CMAKE_SOURCE_DIR}/python/dofulator)
  add_test(
    NAME c_dofulator_in_tree_setup
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:c_dofulator> ${PYTHON_SRC_DIR}
  )
  add_test(
    NAME c_dofulator_in_tree_cleanup
    COMMAND ${CMAKE_COMMAND} -E remove ${PYTHON_SRC_DIR}/$<TARGET_FILE_NAME:c_dofulator>
  )
  set_tests_properties(c_dofulator_in_tree_setup   PROPERTIES FIXTURES_SETUP f_c_dofulator_in_tree)
  set_tests_properties(c_dofulator_in_tree_cleanup PROPERTIES FIXTURES_CLEANUP f_c_dofulator_in_tree)
  set_tests_properties(python_tests PROPERTIES FIXTURES_REQUIRED f_c_dofulator_in_tree)
else()
  message(WARNING "Could not find Pytest! Tests for python interface have been disabled.")
endif()
