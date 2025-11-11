# CMake options=================================================================
cmake_minimum_required(VERSION 3.14) 
# Set your compiler version
set(CMAKE_C_COMPILER gcc-13)
set(CMAKE_CXX_COMPILER g++-13)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
project(GEM)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include(FetchContent)

# Use default FetchContent cache directory (build/_deps)
# set(FETCHCONTENT_BASE_DIR /data/ssalehinasab/cmake_deps)

file(GLOB SRC "src/*.cpp")
file(GLOB CC_SOURCES "thirdparty/plink-2.0/*.cc")
file(GLOB C_SOURCES 
    "thirdparty/zstd-1.5.5/lib/common/*.c"
    "thirdparty/zstd-1.5.5/lib/compress/*.c"
    "thirdparty/zstd-1.5.5/lib/decompress/*.c"
    "thirdparty/libdeflate-1.18/lib/*.c"
    "thirdparty/libdeflate-1.18/lib/x86/*.c"
    "thirdparty/libdeflate-1.18/lib/arm/*.c"
)
# Concatenate the two lists of source files
set(SOURCES ${SRC} ${CC_SOURCES} ${C_SOURCES})
# Fetch fmt ======================================================
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(FMT_POSITION_INDEPENDENT_CODE ON CACHE BOOL "" FORCE) #required for shared libraries forMygen.so
FetchContent_Declare(fmt
 GIT_REPOSITORY https://github.com/fmtlib/fmt.git
 GIT_TAG 11.0.2
)
FetchContent_MakeAvailable(fmt)
#PYBIND11=================================================================
set(PYBIND11_PYTHON_VERSION 3.11)
# set(ENV{LD_LIBRARY_PATH} "$ENV{HOME}/local/Python-3.11.6/")
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.12.0  # or latest stable
)
FetchContent_MakeAvailable(pybind11)


# Fetch Armadillo=================================================================
FetchContent_Declare(
  armadillo
  GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code
  GIT_TAG 14.0.x
)

FetchContent_GetProperties(armadillo)
if(NOT armadillo_POPULATED)
  FetchContent_Populate(armadillo)
  #add_library(armadillo INTERFACE)
  #target_include_directories(armadillo INTERFACE ${armadillo_SOURCE_DIR}/include)
  set(ARMA_USE_OPENMP OFF CACHE BOOL "Disable OpenMP in Armadillo")
  add_definitions(-DARMA_DONT_USE_OPENMP)
  add_subdirectory(${armadillo_SOURCE_DIR} ${armadillo_BINARY_DIR})
endif()
# Fetch Eigen=====================================================
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG 3.4.0
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    add_library(eigen INTERFACE)
    target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
endif()
#Fetch Boost better to install Boost manually=================================================================
# FetchContent_Declare(
#  Boost
#  GIT_REPOSITORY https://github.com/boostorg/boost.git
#  GIT_TAG boost-1.85.0
#  GIT_SHALLOW TRUE
# )
# FetchContent_MakeAvailable(Boost)
# set(BOOST_INCLUDE_LIBRARIES thread filesystem system program_options asio algorithm iostreams math smart_ptr format)
# set(BOOST_ENABLE_CMAKE ON)

# #include(FetchContent)
# FetchContent_Declare(
#   Boost
#   GIT_REPOSITORY https://github.com/boostorg/boost.git
#   GIT_TAG boost-1.79.0
#   GIT_SHALLOW TRUE
# )
# FetchContent_MakeAvailable(Boost)

# Use system Boost installed in Docker
# set(BOOST_ROOT /data/ssalehinasab/boost_1_74_0_install)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.74.0 REQUIRED COMPONENTS program_options thread system filesystem)

include_directories(${Boost_INCLUDE_DIRS})
#find BLAS/LAPACK=================================================================
# message(STATUS "Setting BLAS vendor to Intel10_64lp.")
# set(BLA_VENDOR Intel10_64lp)
# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)
#Set MKL=================================================================
# Use MKL from Docker environment (set by Dockerfile)
if(NOT DEFINED MKLROOT)
  set(MKLROOT $ENV{MKLROOT})
endif()
if(NOT MKLROOT)
  set(MKLROOT "/opt/intel/oneapi/mkl/latest")
endif()
message("MKLROOT: ${MKLROOT}")

# Specify include directories for MKL
include_directories(${MKLROOT}/include)
message(STATUS "Setting include directory for MKL to ${MKLROOT}/include.")

# Specify library directories for MKL
link_directories(${MKLROOT}/lib/intel64)
message(STATUS "Setting library directory for MKL to ${MKLROOT}/lib/intel64.")

# Define EIGEN_USE_MKL_ALL for Eigen library
add_definitions(-DEIGEN_USE_MKL_ALL)
# Define ARMA to use MKL
add_definitions(-DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_DONT_USE_WRAPPER)

message(STATUS "Defining EIGEN_USE_MKL_ALL for Eigen  and Armadillo libraries.")

# Set MKL libraries
set(MKL_LIBRARIES
    ${MKLROOT}/lib/intel64/libmkl_intel_lp64.so
    ${MKLROOT}/lib/intel64/libmkl_sequential.so
    ${MKLROOT}/lib/intel64/libmkl_core.so
    pthread
    m
    dl
)

# Tell CMake that BLAS/LAPACK = MKL
set(BLAS_LIBRARIES ${MKL_LIBRARIES})
set(LAPACK_LIBRARIES ${MKL_LIBRARIES})

message(STATUS "Using MKL for BLAS/LAPACK: ${MKL_LIBRARIES}")


#GMP=================================================================
# GMP and MPFR - use system-installed libraries in Docker
# GMP
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY NAMES gmp libgmp)
if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
  set(GMP_LIBRARIES ${GMP_LIBRARY})
  message(STATUS "Found GMP: ${GMP_LIBRARY}")
else()
  message(WARNING "GMP not found")
endif()

# MPFR
find_path(MPFR_INCLUDE_DIR mpfr.h)
find_library(MPFR_LIBRARY NAMES mpfr libmpfr)
if(MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
  set(MPFR_LIBRARIES ${MPFR_LIBRARY})
  message(STATUS "Found MPFR: ${MPFR_LIBRARY}")
else()
  message(WARNING "MPFR not found")
endif()

if(GMP_INCLUDE_DIR AND MPFR_INCLUDE_DIR)
  include_directories(${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
endif()

message(STATUS "Using GMP include dir: ${GMP_INCLUDE_DIR}")
message(STATUS "Using GMP library: ${GMP_LIBRARIES}")
message(STATUS "Using MPFR include dir: ${MPFR_INCLUDE_DIR}")
message(STATUS "Using MPFR library: ${MPFR_LIBRARIES}")

#SuiteSparse=================================================================
# --- SuiteSparse (use system-installed libraries) -------------------

# Try to find system-installed SuiteSparse
find_library(CHOLMOD_LIB cholmod PATHS /usr/local/lib NO_DEFAULT_PATH)
find_library(CXSPARSE_LIB cxsparse PATHS /usr/local/lib NO_DEFAULT_PATH)
find_library(SPQR_LIB spqr PATHS /usr/local/lib NO_DEFAULT_PATH)
find_library(UMFPACK_LIB umfpack PATHS /usr/local/lib NO_DEFAULT_PATH)
find_library(CAMD_LIB camd PATHS /usr/local/lib NO_DEFAULT_PATH)
find_library(CCOLAMD_LIB ccolamd PATHS /usr/local/lib NO_DEFAULT_PATH)

if(CHOLMOD_LIB AND CXSPARSE_LIB AND SPQR_LIB AND UMFPACK_LIB AND CAMD_LIB AND CCOLAMD_LIB)
    message(STATUS "Using system-installed SuiteSparse libraries")
    set(SUITESPARSE_LIBRARIES
        ${CHOLMOD_LIB}
        ${CXSPARSE_LIB}
        ${SPQR_LIB}
        ${UMFPACK_LIB}
        ${CAMD_LIB}
        ${CCOLAMD_LIB}
        /usr/local/lib/libcolamd.a
        /usr/local/lib/libamd.a
        /usr/local/lib/libsuitesparseconfig.a
    )
    set(SUITESPARSE_INCLUDE_DIR /usr/local/include/suitesparse)
else()
    message(FATAL_ERROR "System-installed SuiteSparse libraries not found")
endif()

# Set executable for Python module ============================================
add_library(Mygen MODULE ${SOURCES})

target_include_directories(Mygen PRIVATE
    ${SUITESPARSE_INCLUDE_DIR}
)
target_include_directories(Mygen PRIVATE ${pybind11_INCLUDE_DIRS})
message(STATUS "pybind11 include dirs: ${pybind11_INCLUDE_DIRS}")
target_link_libraries(Mygen PRIVATE pybind11::module
                  fmt
                  eigen
                  armadillo 
                  Boost::thread
                  Boost::system
                  Boost::filesystem
                  Boost::program_options
                  # nlohmann_json::nlohmann_json
                  ${SUITESPARSE_LIBRARIES}
                  ${MKL_LIBRARIES}
                  m
                  dl)
                  

set_target_properties(Mygen PROPERTIES
                      PREFIX ""  # remove the default 'lib' prefix
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/pymodules"
                     )

# Set compiler optimization flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -funroll-loops -ffast-math  -fno-math-errno -DNDEBUG")
