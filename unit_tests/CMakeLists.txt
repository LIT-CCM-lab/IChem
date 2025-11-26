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
