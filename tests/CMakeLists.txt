# Metatarget for building tests
add_custom_target(tests)

add_subdirectory(unit)
add_subdirectory(integration)

if (NOT DOF_LIBRARY_ONLY)
  add_subdirectory(frontend)
endif()

if (Python_FOUND)
  add_subdirectory(python)
endif()
