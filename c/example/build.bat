rem Define variables
set COMPILE_FLAGS=-O0
set BUILD_DIR=build
set SRC_DIR=..\..\src

rem Create build directory if it does not exist
if not exist %BUILD_DIR% (
    mkdir %BUILD_DIR%
)

rem Compile Fortran files to object files
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack_kinds.F90 -o %BUILD_DIR%\odrpack_kinds.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\lpkbls.f          -o %BUILD_DIR%\lpkbls.o
gfortran %COMPILE_FLAGS% -c %SRC_DIR%\odrpack.f90       -o %BUILD_DIR%\odrpack.o
gfortran %COMPILE_FLAGS% -c ..\odrpack_c.f90            -o %BUILD_DIR%\odrpack_c.o

rem Create the static library
ar rcs %BUILD_DIR%\libodrpack.a %BUILD_DIR%\odrpack_kinds.o %BUILD_DIR%\lpkbls.o %BUILD_DIR%\odrpack.o %BUILD_DIR%\odrpack_c.o

rem Compile C example
gcc %COMPILE_FLAGS% example5.c -L%BUILD_DIR% -lodrpack -lgfortran -o %BUILD_DIR%\example5.exe