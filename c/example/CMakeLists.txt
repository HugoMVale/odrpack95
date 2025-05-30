set(examples 3 5)

get_filename_component(PROJECT_SOURCE_DIR_ABS "${CMAKE_SOURCE_DIR}" ABSOLUTE)

foreach(example IN LISTS examples)
    set(exename "example${example}c")
    set(sourcefile "example${example}.c")

    add_executable(${exename} ${sourcefile})
    target_link_libraries(${exename} PRIVATE odrpack95)
endforeach()
