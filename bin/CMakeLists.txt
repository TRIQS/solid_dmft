set(scripts solid_dmft) 

foreach(script ${scripts})
    configure_file(${script}.py ${script} @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endforeach()
