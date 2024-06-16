#!/bin/bash

# Define variables
COMPILE_FLAGS="-O0"
BUILD_DIR="build"
SRC_DIR="../../src"

# Create build directory if it does not exist
if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
fi

# Compile Fortran files to object files
gfortran $COMPILE_FLAGS -c $SRC_DIR/odrpack_kinds.F90 -o $BUILD_DIR/odrpack_kinds.o
gfortran $COMPILE_FLAGS -c $SRC_DIR/lpkbls.f          -o $BUILD_DIR/lpkbls.o
gfortran $COMPILE_FLAGS -c $SRC_DIR/odrpack.f90       -o $BUILD_DIR/odrpack.o
gfortran $COMPILE_FLAGS -c ../odrpack_c.f90           -o $BUILD_DIR/odrpack_c.o

# Create the static library
ar rcs $BUILD_DIR/libodrpack.a $BUILD_DIR/odrpack_kinds.o $BUILD_DIR/lpkbls.o $BUILD_DIR/odrpack.o $BUILD_DIR/odrpack_c.o

# Compile C example
gcc $COMPILE_FLAGS example5.c -L$BUILD_DIR -lodrpack -lgfortran -o $BUILD_DIR/example5
