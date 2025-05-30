cmake_minimum_required(VERSION 3.16)
project(odrpack95 VERSION 1.0.1 LANGUAGES C Fortran)

# Enable testing
enable_testing()

# Build type and Fortran standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_Fortran_STANDARD 90)

# Option flags for examples and tests
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTS "Build tests" OFF)

# Include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/c/include)

# BLAS detection
set(BLAS_FOUND FALSE)
if(APPLE)
    message(STATUS "Using Accelerate framework as BLAS library.")
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    if(ACCELERATE_FRAMEWORK)
        set(BLAS_LIBRARIES ${ACCELERATE_FRAMEWORK})
        set(BLAS_FOUND TRUE)
    endif()
else()
    find_package(PkgConfig)
    pkg_check_modules(OPENBLAS openblas)
    if(OPENBLAS_FOUND)
        set(BLAS_LIBRARIES ${OPENBLAS_LIBRARIES})
        set(BLAS_FOUND TRUE)
    else()
        find_library(OPENBLAS_LIB openblas)
        if(OPENBLAS_LIB)
            set(BLAS_LIBRARIES ${OPENBLAS_LIB})
            set(BLAS_FOUND TRUE)
        endif()
    endif()
endif()

# If no BLAS found, fall back to reference implementation
if(NOT BLAS_FOUND)
    message(STATUS "No specific BLAS library detected. Falling back to reference BLAS.")
    set(BLAS_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/extern/blas.f)
else()
    set(BLAS_SOURCE "")
endif()

# ODRPACK source files
set(ODRPACK_SRC
    src/linpack.f
    src/blas_interfaces.f90
    src/odrpack.f90
    src/odrpack_kinds.F90
    src/odrpack_core.f90
    src/odrpack_reports.f90
    src/odrpack_capi.f90
    ${BLAS_SOURCE}
)

# Create static library
add_library(odrpack95 STATIC ${ODRPACK_SRC})
target_link_libraries(odrpack95 PRIVATE ${BLAS_LIBRARIES})

# Expose the module directory to targets that link with odrpack95
target_include_directories(odrpack95 PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Build examples if requested
if(BUILD_EXAMPLES)
    add_subdirectory(example)
    add_subdirectory(c/example)
endif()

# Build tests if requested
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
