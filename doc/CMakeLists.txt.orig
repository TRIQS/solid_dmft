

# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)


# ---------------------------------
# Top Sphinx target
# ---------------------------------
if(NOT DEFINED SPHINXBUILD_EXECUTABLE)
  find_package(Sphinx)
endif()


# Sphinx has internal caching, always run it
add_custom_target(${PROJECT_NAME}_docs_sphinx ALL)
add_custom_command(
  TARGET ${PROJECT_NAME}_docs_sphinx
  COMMAND cd ./input_output/DMFT_input && bash generate_doc_from_comments.sh && cd ../..
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:${PROJECT_BINARY_DIR}/python/solid_dmft:$ENV{PYTHONPATH} ${SPHINXBUILD_EXECUTABLE} -c . -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/${PROJECT_NAME}
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|js|xsl|css|py|txt|inv|bib)$"
  PATTERN "_*"
)
