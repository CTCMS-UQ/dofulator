
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# Cythonize wrapper into .c file
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/c_dofulator.c
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/c_dofulator.pyx
    ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/c_dofulator.pxd
  VERBATIM
  COMMAND
    Python::Interpreter -m cython ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/c_dofulator.pyx
    --output-file ${CMAKE_CURRENT_BINARY_DIR}/c_dofulator.c
)

# Compile generated .c file
# Static linkage of libdofulator here gets blas libraries linked into the final packaged c_dofulator shared lib
python_add_library(c_dofulator MODULE ${CMAKE_CURRENT_BINARY_DIR}/c_dofulator.c WITH_SOABI)
set_target_properties(c_dofulator PROPERTIES EXCLUDE_FROM_ALL YES)
target_include_directories(c_dofulator PRIVATE ${Python_NumPy_INCLUDE_DIRS} ${SRC_DIR})
target_compile_definitions(c_dofulator PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_link_libraries(c_dofulator PRIVATE dofulator)

file(GLOB PY_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/dofulator/*.*)
install(TARGETS c_dofulator DESTINATION dofulator/)
install(FILES ${PY_SRCS} DESTINATION dofulator/)
