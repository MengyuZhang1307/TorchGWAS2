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
set(ENV{LD_LIBRARY_PATH} "$ENV{HOME}/local/Python-3.11.6/")
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
#Fetch Boost=================================================================
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

#find Boost=================================================================
find_package(Boost 1.74.0 REQUIRED COMPONENTS program_options thread system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
#find BLAS/LAPACK=================================================================
# message(STATUS "Setting BLAS vendor to Intel10_64lp.")
# set(BLA_VENDOR Intel10_64lp)
# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)
#Set MKL=================================================================
set(MKLROOT "/opt/intel/compilers_and_libraries_2020.2.254/linux/mkl")
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
    mkl_intel_lp64       
    mkl_sequential    
    mkl_core  
    mkl_def
    mkl_vml_def
    pthread
    m
    dl
)


#SuiteSparse=================================================================
#If SuiteSparse is installed in default directory uncomment below line and comment all the lines related to the SuiteSparse
#find_package(SuiteSparse REQUIRED)
# Define SuiteSparse installation directory
# Include the SuiteSparse headers compiled with suitesparse_7.6.0_no_omp replace your path
include_directories("/home/ssalehinasab/local/suitesparse_7.6.0_no_omp/include/suitesparse/")

# Link the SuiteSparse libraries
link_directories("/home/ssalehinasab/local/suitesparse_7.6.0_no_omp/lib64")

# Add the SuiteSparse libraries to the target
set(SUITESPARSE_LIBRARIES
    libcholmod.a 
    libcxsparse.a
    libspqr.a 
    libumfpack.a 
    libcamd.a 
    libccolamd.a 
    libcolamd.a 
    libamd.a 
    libsuitesparseconfig.a 
)

# Set executable for Python module ================================
set(SOURCES ${SRC} ${CC_SOURCES} ${C_SOURCES})
add_library(Mygen MODULE ${SOURCES})
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
                  ${SUITESPARSE_LIBRARIES}
                  ${MKL_LIBRARIES})
set_target_properties(Mygen PROPERTIES PREFIX "")

# Set compiler optimization flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -funroll-loops -ffast-math  -fno-math-errno -DNDEBUG")
