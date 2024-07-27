#!/bin/bash

# Define variables
COMPILE_FLAGS="-O0"
BUILD_DIR="build"
SRC_DIR="../../src"

# Create build directory if it does not exist
if [ ! -d "$BUILD_DIR" ]; then
    mkdir "$BUILD_DIR"
fi

# Compile Fortran files to object files
gfortran $COMPILE_FLAGS -c "$SRC_DIR/odrpack_kinds.F90"   -o "$BUILD_DIR/odrpack_kinds.o"
gfortran $COMPILE_FLAGS -c "$SRC_DIR/blas.f"              -o "$BUILD_DIR/blas.o"
gfortran $COMPILE_FLAGS -c "$SRC_DIR/linpack.f"           -o "$BUILD_DIR/linpack.o"
gfortran $COMPILE_FLAGS -c "$SRC_DIR/odrpack_core.f90"    -o "$BUILD_DIR/odrpack_core.o"
gfortran $COMPILE_FLAGS -c "$SRC_DIR/odrpack_reports.f90" -o "$BUILD_DIR/odrpack_reports.o"
gfortran $COMPILE_FLAGS -c "$SRC_DIR/odrpack.f90"         -o "$BUILD_DIR/odrpack.o"
gfortran $COMPILE_FLAGS -c "$SRC_DIR/odrpack_capi.f90"    -o "$BUILD_DIR/odrpack_capi.o"

# Create the static library
ar rcs "$BUILD_DIR/libodrpack.a" "$BUILD_DIR/odrpack_kinds.o" "$BUILD_DIR/blas.o" "$BUILD_DIR/linpack.o" "$BUILD_DIR/odrpack_core.o" "$BUILD_DIR/odrpack_reports.o" "$BUILD_DIR/odrpack.o" "$BUILD_DIR/odrpack_capi.o"

# Compile C example
gcc $COMPILE_FLAGS example5.c -L"$BUILD_DIR" -lodrpack -lgfortran -o "$BUILD_DIR/example5.exe"
