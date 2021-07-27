set(scripts solid_dmft symmetrize_gamma_file update_solid_config update_solid_h5)

foreach(script ${scripts})
    configure_file(${script}.py ${script} @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endforeach()
