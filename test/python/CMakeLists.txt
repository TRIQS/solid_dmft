# all pytest unittests
set (all_pytests
     test_afm_mapping
     test_interaction_hamiltonian
     test_manipulate_chemical_potential.py
     test_observables.py
     test_read_config.py
     test_update_dmft_config.py
     test_update_results_h5.py
    )

foreach(test ${all_pytests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)

  add_test(NAME ${test_name}
           COMMAND ${TRIQS_PYTHON_EXECUTABLE} -m pytest -vv ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})

  set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH})
endforeach()

# ------------------------------#

# all other tests
set(all_tests
    test_convergence
    test_matheval
    test_plot_correlated_bands
    test_respack_sfo
    test_gw_embedding
    )

  # copy input files for PCB test from documentation tutorial
  set (test_dir ${CMAKE_CURRENT_BINARY_DIR})
  foreach(file svo.wout svo_example.h5 svo_hr.dat)
    FILE(COPY ../../doc/tutorials/correlated_bandstructure/${file} DESTINATION ${test_dir})
  endforeach()

  # copy reference data for PCB test
  FILE(COPY test_pcb_ref.h5 DESTINATION ${test_dir})

  # copy reference data for respack test
  FILE(COPY respack_sfo_data DESTINATION ${test_dir})

  # copy reference data for respack test
  FILE(COPY gw_embedding DESTINATION ${test_dir})
  
foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME ${test_name} COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH})
endforeach()

# ------------------------------#

# integration tests
set (integration_tests
     svo_hubbardI_basic
     svo_hartree
     lno_hubbardI_mag
     svo_cthyb_basic_gl
     svo_cthyb_basic_tf
     svo_ctseg_dyn
     )

# copy other test files
FILE(COPY UIJKL DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(test ${integration_tests})
  set (test_dir ${CMAKE_CURRENT_BINARY_DIR}/${test})

  foreach(file dmft_config.ini inp.h5 ref.h5 test.py)
    FILE(COPY ${test}/${file} DESTINATION ${test_dir})
  endforeach()

  add_test(NAME ${test}
           #COMMAND bash ${test}.sh
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${TRIQS_PYTHON_EXECUTABLE} test.py
           WORKING_DIRECTORY ${test_dir}
          )

  set_property(TEST ${test} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH})
endforeach()

# ------------------------------#

# MaxEnt test on svo_hubbardI_basic
set (test svo_hubbardI_basic)

foreach(file test_maxent.py helper.py)
    FILE(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

add_test(NAME maxent_${test}
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${TRIQS_PYTHON_EXECUTABLE} test_maxent.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )

set_property(TEST maxent_${test} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH})
