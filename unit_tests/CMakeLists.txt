# ICTools/ifp_test.cpp test
add_executable(ifp_unit_tests
    ICTools/ifp_test.cpp
)

target_link_libraries(ifp_unit_tests
    PRIVATE
        ichem_core
)

add_test(
    NAME ICTools_ifp_unit
    COMMAND ifp_unit_tests
)

# ICMole/molecule_test.cpp: molecule editing on test/dataset ligands
add_executable(molecule_unit_tests
    ICMole/molecule_test.cpp
)

target_link_libraries(molecule_unit_tests
    PRIVATE
        ichem_core
)

target_compile_definitions(molecule_unit_tests
    PRIVATE
        ICHEM_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

add_test(
    NAME ICMole_molecule_unit
    COMMAND molecule_unit_tests
)
