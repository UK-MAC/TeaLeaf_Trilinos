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
rm -rf CMakeCache.txt CMakeFiles

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
# Added -g flag to enable debugging/profiling
#
# Optional MKL BLAS / LAPACK setting   
#    -D TPL_ENABLE_MKL:BOOL=ON \
#    -D MKL_LIBRARY_DIRS:STRING="/opt/intel/Compiler/16.0/current/mkl/lib/intel64" \
#    -D BLAS_LIBRARY_DIRS:STRING="/opt/intel/Compiler/16.0/current/mkl/lib/intel64;/opt/intel/Compiler/16.0/current/compiler/lib/intel64;/usr/lib64" \
#    -D BLAS_LIBRARY_NAMES:STRING="mkl_intel_lp64; mkl_intel_thread; mkl_core; iomp5; pthread" \
#    -D LAPACK_LIBRARY_DIRS:STRING="/opt/intel/Compiler/16.0/current/mkl/lib/intel64;/opt/intel/Compiler/16.0/current/compiler/lib/intel64;/usr/lib64" \
#    -D LAPACK_LIBRARY_NAMES:STRING="mkl_intel_lp64; mkl_intel_thread; mkl_core; iomp5; pthread" \
#   
# export PATH=/opt/gcc/gcc-5.2.0/bin/:$PATH
# export LD_LIBRARY_PATH=/opt/gcc/gcc-5.2.0/lib64:$LD_LIBRARY_PATH
# Note LAPACK 3.6+ is not supported yet

## Build option for Trilinos 13.0.1 using gcc/7.3.1 + spectrum-mpi-2020.08.19 with cuda/10.1.243 for V100 GPUs

TRILINOS_PATH=<Trilinos Source DIR>
LOCALBASE=<MPI INSTALL DIR>
PREFIX=~/TeaLeaf_Trilinos/libs/trilinos
MPICH_CXX=`readlink -f ../packages/kokkos/config/nvcc_wrapper`
#export NVCC_WRAPPER_DEFAULT_COMPILER=g++
#export CUDA_LAUNCH_BLOCKING=1


cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D MPI_BASE_DIR="${LOCALBASE}" \
    -D MPI_CXX_COMPILER:FILEPATH="${LOCALBASE}/bin/mpic++"    \
    -D CMAKE_CXX_COMPILER:FILEPATH="${TRILINOS_PATH}/packages/kokkos/bin/nvcc_wrapper" \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
    -D Trilinos_ENABLE_MueLu:BOOL=ON \
    -D Trilinos_ENABLE_ML:BOOL=ON \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Trilinos_ENABLE_Xpetra:BOOL=ON \
    -DTpetra_ENABLE_DEPRECATED_CODE:BOOL=ON \
    -D Trilinos_ENABLE_OpenMP:BOOL=ON \
    -D Tpetra_INST_SERIAL:BOOL=ON \
    -D Tpetra_INST_OPENMP:BOOL=ON \
    -D Tpetra_INST_INT_INT:BOOL=ON \
    -DXpetra_ENABLE_DEPRECATED_CODE:BOOL=ON \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
    -D Trilinos_ENABLE_TESTS:BOOL=OFF \
    -D MueLu_ENABLE_EXAMPLES:BOOL=OFF \
    -D MueLu_ENABLE_TESTS:BOOL=OFF \
    -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
    -D Trilinos_ENABLE_COMPLEX:BOOL=OFF \
    -D CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON \
    -D Trilinos_ENABLE_CXX11:BOOL=ON \
    -D Trilinos_CXX11_FLAGS="--expt-extended-lambda -std=c++11" \
    -D Trilinos_ENABLE_Fortran:BOOL=OFF\
    -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
    -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D TPL_ENABLE_CUDA:BOOL=ON \
    -D TPL_ENABLE_CUSPARSE:BOOL=ON \
    -D TPL_ENABLE_CUSPARSE:BOOL=ON \
    -D TPL_ENABLE_HWLOC:BOOL=OFF \
    -D TPL_ENABLE_BLAS:BOOL=ON \
    -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D Kokkos_ENABLE_OPENMP:BOOL=ON \
    -D Kokkos_ENABLE_CUDA_UVM:BOOL=ON \
    -D Kokkos_ENABLE_CUDA:BOOL=ON \
    -D Kokkos_ENABLE_CUDA_LAMBDA:BOOL=ON \
    -D Kokkos_ARCH_POWER9=ON \
    -D Kokkos_ARCH_VOLTA70=ON \
    $EXTRA_ARGS \
    ../


make -j32
make install -j32
