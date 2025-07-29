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


set(SRC src/ReadFiles.cpp
        src/ReadBGEN.cpp
        src/main.cpp
        src/RunPipeline.cpp)

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
#=================================================================
#=================================================================

# Download and link pybind11############################################################
#If it is not able to find right python then do:
# on terminal set:
# unset PYTHONHOME
# unset PYTHONPATH
#Can run camke as below:
#cmake -B build -S . -DPython3_EXECUTABLE=$HOME/local/python3.11/bin/python3.11
#or set variable inside cmake as:set(Python3_EXECUTABLE "$ENV{HOME}/local/python3.11/bin/python3.11")
#Since pybind cannot find the right python version so we do:
# If you do not have python issue comment the following 2 lines
set(PYBIND11_PYTHON_VERSION 3.11)
set(Python3_EXECUTABLE "$ENV{HOME}/local/python3.11/bin/python3.11")
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.12.0  # or latest stable
)
FetchContent_MakeAvailable(pybind11)
# Set executable for Python module ================================
set(SRC_PY src/ReadFiles.cpp
                src/ReadBGEN.cpp
                src/PyInterface.cpp
                src/RunPipeline.cpp
                ${CC_SOURCES} ${C_SOURCES})
add_library(Mygen MODULE ${SRC_PY})
target_link_libraries(Mygen PRIVATE pybind11::module)
set_target_properties(Mygen PROPERTIES PREFIX "")

# Add the executable with the concatenated source files for C++ binary
# add_executable(GEM ${SOURCES} )
# target_link_libraries(GEM PRIVATE pthread)

# Set compiler optimization flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
