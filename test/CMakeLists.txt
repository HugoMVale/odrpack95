set(tests
    error_detection
    solution
)

get_filename_component(PROJECT_SOURCE_DIR_ABS "${CMAKE_SOURCE_DIR}" ABSOLUTE)

foreach(t IN LISTS tests)
    set(testname "test-${t}")
    set(sourcefile "test_${t}.f90")

    add_executable(${testname} ${sourcefile})
    target_link_libraries(${testname} PRIVATE odrpack95)

    add_test(NAME ${t} COMMAND ${testname})

    set_tests_properties(${t} PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR_ABS}"
    )
endforeach()
