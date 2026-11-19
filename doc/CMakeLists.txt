# Fail early if a doc dependency is missing (import names, not distribution names)
find_package(Python COMPONENTS Interpreter REQUIRED)

foreach(_module numpydoc nbsphinx myst_parser linkify_it sphinx_rtd_theme matplotlib sphinx_autorun)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import ${_module}"
    RESULT_VARIABLE _rc OUTPUT_QUIET ERROR_QUIET)
  if(_rc)
    message(FATAL_ERROR "Building the documentation requires the python module '${_module}' (see requirements.txt).")
  endif()
endforeach()

# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
if(NOT DEFINED SPHINXBUILD_EXECUTABLE)
  find_package(Sphinx)
endif()

# Sphinx has internal caching, always run it
add_custom_target(${PROJECT_NAME}_docs_sphinx ALL)
# rst pages of the input parameters, generated from io_tools/documentation.txt
add_custom_command(
  TARGET ${PROJECT_NAME}_docs_sphinx POST_BUILD
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${TRIQS_PYTHON_EXECUTABLE} -m solid_dmft.io_tools.documentation ${CMAKE_CURRENT_SOURCE_DIR}/input_output/DMFT_input
)
add_custom_command(
  TARGET ${PROJECT_NAME}_docs_sphinx POST_BUILD
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:${PROJECT_BINARY_DIR}/python/solid_dmft:$ENV{PYTHONPATH} ${SPHINXBUILD_EXECUTABLE} -c . -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/${PROJECT_NAME}
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|ico|js|xsl|css|py|txt|inv|bib|ttf|woff2|eot|sh)$"
  PATTERN "_*"
)
