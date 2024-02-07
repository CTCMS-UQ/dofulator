add_subdirectory(unit)

if (NOT DOF_LIBRARY_ONLY)
  add_subdirectory(integration)
endif()
