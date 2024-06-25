# Output from running the frontend on each .xyz file stored as a .dof file
file(GLOB FRONTEND_TEST_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/*.xyz)

add_custom_target(frontend_test_adopt_all)

foreach(XYZ IN ITEMS ${FRONTEND_TEST_INPUT})
  get_filename_component(FRONTEND_TEST ${XYZ} NAME_WLE)
  set(FRONTEND_TEST_REF ${CMAKE_CURRENT_SOURCE_DIR}/${FRONTEND_TEST}.dof)
  set(FRONTEND_TEST_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${FRONTEND_TEST}.dof)
  # TODO: Make this platform independent
  add_test(
    NAME frontend::${FRONTEND_TEST}
    COMMAND bash -c "${CMAKE_BINARY_DIR}/dof ${XYZ} -o ${FRONTEND_TEST_OUTPUT} && diff ${FRONTEND_TEST_REF} ${FRONTEND_TEST_OUTPUT}"
  )

  # Target to write .dof file using current output
  add_custom_target(
    frontend_test_adopt_${FRONTEND_TEST}
    DEPENDS dof ${XYZ}
    USES_TERMINAL
    COMMAND dof ${XYZ} -o ${FRONTEND_TEST_REF}
  )
  add_dependencies(frontend_test_adopt_all frontend_test_adopt_${FRONTEND_TEST})
endforeach()
