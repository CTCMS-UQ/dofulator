file(GLOB INT_TEST_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/*.xyz)

add_custom_target(inttest_adopt_all)

foreach(XYZ IN ITEMS ${INT_TEST_INPUT})
  get_filename_component(INT_TEST ${XYZ} NAME_WLE)
  set(INT_TEST_REF ${CMAKE_CURRENT_SOURCE_DIR}/${INT_TEST}.dof)
  set(INT_TEST_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${INT_TEST}.dof)
  # TODO: Make this platform independent
  add_test(
    NAME ${XYZ}
    COMMAND bash -c "${CMAKE_BINARY_DIR}/dof ${XYZ} -o ${INT_TEST_OUTPUT} && diff ${INT_TEST_REF} ${INT_TEST_OUTPUT}"
  )

  add_custom_target(
    inttest_adopt_${INT_TEST}
    DEPENDS dof ${XYZ}
    USES_TERMINAL
    COMMAND dof ${XYZ} -o ${INT_TEST_REF}
  )
  add_dependencies(inttest_adopt_all inttest_adopt_${INT_TEST})
endforeach()
