
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# Cythonize wrapper into .c file
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/c_dofulator.c
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/c_dofulator.pyx
    ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/c_dofulator.pxd
  VERBATIM
  COMMAND
    Python::Interpreter -m cython -3 ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/c_dofulator.pyx
    --output-file ${CMAKE_CURRENT_BINARY_DIR}/c_dofulator.c
)

# Compile generated .c file
# Static linkage of libdofulator here gets blas libraries linked into the final packaged c_dofulator shared lib
python_add_library(c_dofulator MODULE ${CMAKE_CURRENT_BINARY_DIR}/c_dofulator.c WITH_SOABI)
set_target_properties(c_dofulator PROPERTIES EXCLUDE_FROM_ALL YES)
target_include_directories(c_dofulator PRIVATE ${Python_NumPy_INCLUDE_DIRS} ${SRC_DIR})
target_compile_definitions(c_dofulator PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_link_libraries(c_dofulator PRIVATE dofulator)

# FindPython sometimes gives incorrect SOABI suffix. Use SKBUILD to correct it if available.
if (SKBUILD)
  if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set_property (TARGET ${name} PROPERTY SUFFIX ".${SKBUILD_SOABI}.pyd")
  else()
    set_property (TARGET ${name} PROPERTY SUFFIX ".${SKBUILD_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
  endif()
endif()

# Avoid double extension bug
get_target_property(SOABI_SUFFIX c_dofulator SUFFIX)
get_filename_component(CDOF_EXT1 ${SOABI_SUFFIX} LAST_EXT)
get_filename_component(CDOF_NOEXT1 ${SOABI_SUFFIX} NAME_WLE)
get_filename_component(CDOF_EXT2 ${CDOF_NOEXT1} LAST_EXT)
if (CDOF_EXT1 STREQUAL CDOF_EXT2)
  set_target_properties(c_dofulator PROPERTIES SUFFIX ${CDOF_NOEXT1})
endif()

# Mark as windows build if needed
if (WIN32)
  target_compile_definitions(c_dofulator PRIVATE MS_WIN32 MS_WINDOWS)
cmake_host_system_information(RESULT IS_64BIT QUERY IS_64BIT)
  if (IS_64BIT)
    target_compile_definitions(c_dofulator PRIVATE MS_WIN64)
  endif()
endif()

# Install details
file(GLOB PY_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/*.*)
install(TARGETS c_dofulator DESTINATION dofulator/)
install(FILES ${PY_SRCS} DESTINATION dofulator/)
