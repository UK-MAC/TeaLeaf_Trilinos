#!/bin/bash

# Youu can invoke this shell script with additional command-line
# arguments.  They will be passed directly to CMake.
#
EXTRA_ARGS=$@

#
# Each invocation of CMake caches the values of build options in a
# CMakeCache.txt file.  If you run CMake again without deleting the
# CMakeCache.txt file, CMake won't notice any build options that have
# changed, because it found their original values in the cache file.
# Deleting the CMakeCache.txt file before invoking CMake will insure
# that CMake learns about any build options you may have changed.
#
rm -f CMakeCache.txt

#
# A sampling of CMake build options.
#
# CMAKE_INSTALL_PREFIX: Where to install Trilinos.
# CMAKE_BUILD_TYPE: DEBUG or RELEASE.
# CMAKE_CXX_COMPILER: The C++ compiler to use when building Trilinos.
# CMAKE_C_COMPILER: The C compiler to use when building Trilinos.
#   Some parts of Trilinos are implemented in C.
# CMAKE_Fortran_COMPILER: The Fortran compiler to use when building
#   Trilinos.  Some parts of Trilinos are implemented in Fortran.
# HAVE_GCC_ABI_DEMANGLE: Setting this option to ON improves 
#   debugging messages.
# CMAKE_VERBOSE_MAKEFILE: Set to OFF (or FALSE) if you prefer a quiet
#   build.
# Trilinos_ENABLE_ALL_PACKAGES: If you like, you can build _all_ of
#   Trilinos, but you don't have to.
# Trilinos_ENABLE_Epetra: If ON, build Epetra.
# Trilinos_ENABLE_Triutils: If ON, bulid Triutils.
# Trilinos_ENABLE_TESTS: If ON, build the tests for all packages that
#   are to be built.
# Trilinos_ENABLE_EXAMPLES: If ON, build the examples for all
#   packages that are to be built.
#

PREFIX=/home/wpg/TeaLeaf/tealeaf_trilinos/libs/trilinos

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D CMAKE_CXX_COMPILER:FILEPATH=mpiicpc \
    -D CMAKE_C_COMPILER:FILEPATH=mpiicc \
    -D CMAKE_Fortran_COMPILER:FILEPATH=mpiifort  \
    -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE \
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
    -D Trilinos_ENABLE_TESTS:BOOL=ON \
    -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D BLAS_LIBRARY_DIRS:FILEPATH=/opt/lapack/3.4.2/intel-13.1.1.163/lib \
    -D LAPACK_LIBRARY_DIRS:FILEPATH=/opt/lapack/3.4.2/intel-13.1.1.163/lib \
    -D CMAKE_CXX_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK" \
    $EXTRA_ARGS \
    ../

make
make install