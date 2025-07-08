cmake_minimum_required(VERSION 3.14) 
# Set your compiler version
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
project(GEM)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include(FetchContent)

# CMake options=================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
# Add the executable with the concatenated source files
add_executable(GEM ${SOURCES} )
# Download and link pybind11############################################################
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.12.0  # or latest stable
)
FetchContent_MakeAvailable(pybind11)
# Set executable for Python module ================================
set(SRC_PY src/ReadFiles.cpp
        src/PyInterface.cpp
        src/main.cpp
        src/RunPipeline.cpp
        ${CC_SOURCES} ${C_SOURCES})
add_library(Mygen MODULE ${SRC_PY})
target_link_libraries(Mygen PRIVATE pybind11::module)
# Set compiler optimization flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")