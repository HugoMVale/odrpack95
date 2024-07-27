rem Define variables
set COMPILE_FLAGS=-O0
set BUILD_DIR=build
set SRC_DIR=..\..\src

rem Create build directory if it does not exist
if not exist %BUILD_DIR% (
    mkdir %BUILD_DIR%
)

rem Compile Fortran files to object files
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack_kinds.F90   -o %BUILD_DIR%\odrpack_kinds.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\blas.f              -o %BUILD_DIR%\blas.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\linpack.f           -o %BUILD_DIR%\linpack.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack_core.f90    -o %BUILD_DIR%\odrpack_core.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack_reports.f90 -o %BUILD_DIR%\odrpack_reports.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack.f90         -o %BUILD_DIR%\odrpack.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack_capi.f90    -o %BUILD_DIR%\odrpack_capi.o

rem Create the static library
ar rcs %BUILD_DIR%\libodrpack.a %BUILD_DIR%\odrpack_kinds.o %BUILD_DIR%\blas.o %BUILD_DIR%\linpack.o %BUILD_DIR%\odrpack_core.o %BUILD_DIR%\odrpack_reports.o %BUILD_DIR%\odrpack.o %BUILD_DIR%\odrpack_capi.o

rem Compile C example
gcc %COMPILE_FLAGS% example5.c -L%BUILD_DIR% -lodrpack -lgfortran -o %BUILD_DIR%\example5.exe