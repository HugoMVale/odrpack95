cmake_minimum_required(VERSION 3.16)
project(odrpack95 VERSION 3.1.0 LANGUAGES C Fortran)

# Enable testing
enable_testing()

# Build type and Fortran standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_Fortran_STANDARD 90)

# Option flags for examples and tests
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTS "Build tests" OFF)
option(BUILD_SHARED "Build a shared library" OFF)
set(BLAS_OPTION "auto" CACHE STRING "BLAS library to use (auto, openblas, or refblas)")
set_property(CACHE BLAS_OPTION PROPERTY STRINGS "auto" "openblas" "refblas")

# Include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/c/include)

# BLAS detection
set(BLAS_FOUND FALSE)
set(BLAS_SOURCE "")

if(BLAS_OPTION STREQUAL "refblas")
    message(STATUS "`Reference BLAS` selected.")
    set(BLAS_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/extern/blas.f)

elseif(BLAS_OPTION STREQUAL "openblas")
    message(STATUS "`OpenBLAS` selected.")
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(OPENBLAS openblas)
    endif()
    
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
    
    if(NOT BLAS_FOUND)
        message(FATAL_ERROR "OpenBLAS library not found. Install OpenBLAS or select another BLAS option.")
    endif()

elseif(BLAS_OPTION STREQUAL "auto")
    message(STATUS "`Auto` BLAS selection enabled.")
    
    if(APPLE)
        find_library(ACCELERATE_FRAMEWORK Accelerate)
        if(ACCELERATE_FRAMEWORK)
            message(STATUS "Using Apple Accelerate framework (BLAS/LAPACK).")
            set(BLAS_LIBRARIES ${ACCELERATE_FRAMEWORK})
            set(BLAS_FOUND TRUE)
        endif()
    endif()
    
    if(NOT BLAS_FOUND)
        foreach(lib IN ITEMS blas openblas)
            find_package(PkgConfig QUIET)
            if(PkgConfig_FOUND)
                pkg_check_modules(BLAS_LIB ${lib})
            endif()
            
            if(BLAS_LIB_FOUND)
                message(STATUS "Detected BLAS library: ${lib}")
                set(BLAS_LIBRARIES ${BLAS_LIB_LIBRARIES})
                set(BLAS_FOUND TRUE)
                break()
            else()
                find_library(BLAS_LIB_PATH ${lib})
                if(BLAS_LIB_PATH)
                    message(STATUS "Detected BLAS library: ${lib}")
                    set(BLAS_LIBRARIES ${BLAS_LIB_PATH})
                    set(BLAS_FOUND TRUE)
                    break()
                endif()
            endif()
        endforeach()
    endif()
    
    if(NOT BLAS_FOUND)
        message(WARNING "No BLAS library detected. Falling back to reference BLAS.")
        set(BLAS_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/extern/blas.f)
    endif()

else()
    message(FATAL_ERROR "Invalid BLAS option: ${BLAS_OPTION}")
endif()

# ODRPACK source files
set(ODRPACK_SRC
    src/linpack.f
    src/blas_interfaces.f90
    src/odrpack.f90
    src/odrpack_kinds.f90
    src/odrpack_core.f90
    src/odrpack_reports.f90
    src/odrpack_capi.f90
    ${BLAS_SOURCE}
)

# Create static or shared library
if(BUILD_SHARED)
   add_library(odrpack95 SHARED ${ODRPACK_SRC})
else()   
   add_library(odrpack95 STATIC ${ODRPACK_SRC})
endif()

# Add GFORTRAN compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "Adding GNU Fortran-specific flags to odrpack95.")
    target_compile_options(odrpack95 PUBLIC 
        -Wno-compare-reals
        -Wno-unused-dummy-argument
    )
endif()

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

# Install the library
include(GNUInstallDirs)
install(TARGETS odrpack95
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
