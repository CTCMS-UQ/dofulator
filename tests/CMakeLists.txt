add_subdirectory(unit)
add_subdirectory(integration)

if (NOT DOF_LIBRARY_ONLY)
  add_subdirectory(frontend)
endif()
