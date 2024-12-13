
# Set up virtual environment with test dependencies
add_custom_command(
  OUTPUT venv
  COMMAND Python::Interpreter -m venv venv
)
set(VENV_Python ${CMAKE_CURRENT_BINARY_DIR}/venv/bin/python)
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/venv.stamp
  DEPENDS venv ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
  COMMAND ${VENV_Python} -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt --upgrade
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/venv.stamp
)
add_custom_target(
  pytest_venv
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/venv.stamp
)

# Use pytest from venv
add_test(
  NAME python_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND ${VENV_Python} -m pytest
)

# Set up fixture to copy c_dofulator library into source tree/clean it up afterwards for pytest
add_dependencies(tests c_dofulator pytest_venv)
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
