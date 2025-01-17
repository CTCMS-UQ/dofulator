# Output from running the frontend on each .xyz file stored as a .dof file
file(GLOB FRONTEND_TEST_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/*.xyz)

add_custom_target(frontend_test_adopt_all)

foreach(XYZ IN ITEMS ${FRONTEND_TEST_INPUT})
  get_filename_component(FRONTEND_TEST ${XYZ} NAME_WLE)
  set(FRONTEND_TEST_REF ${CMAKE_CURRENT_SOURCE_DIR}/${FRONTEND_TEST}.dof)
  set(FRONTEND_TEST_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${FRONTEND_TEST}.dof)
  add_test(
    NAME frontend::run_${FRONTEND_TEST}
    COMMAND $<TARGET_FILE:dof> ${XYZ} -o ${FRONTEND_TEST_OUTPUT}
  )
  add_test(
    NAME frontend::check_${FRONTEND_TEST}
    COMMAND ${CMAKE_COMMAND} -E compare_files ${FRONTEND_TEST_REF} ${FRONTEND_TEST_OUTPUT}
  )
  set_tests_properties(frontend::run_${FRONTEND_TEST}   PROPERTIES FIXTURES_SETUP    f_${FRONTEND_TEST})
  set_tests_properties(frontend::check_${FRONTEND_TEST} PROPERTIES FIXTURES_REQUIRED f_${FRONTEND_TEST})
  # Need to include scipy-openblas32 lib directory in PATH for tests to find the shared object
  if (DOF_USE_SCIPY_OPENBLAS32)
    set_tests_properties(frontend::run_${FRONTEND_TEST} PROPERTIES ENVIRONMENT_MODIFICATION PATH=path_list_prepend:${SCIPY_OPENBLAS_LIB})
  endif()
  add_dependencies(tests dof)

  # Target to write .dof file using current output
  add_custom_target(
    frontend_test_adopt_${FRONTEND_TEST}
    DEPENDS dof ${XYZ}
    USES_TERMINAL
    COMMAND dof ${XYZ} -o ${FRONTEND_TEST_REF}
  )
  add_dependencies(frontend_test_adopt_all frontend_test_adopt_${FRONTEND_TEST})
endforeach()
