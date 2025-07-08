# cmake_minimum_required(VERSION 3.25)
# project(ReadBgen)
# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

# # Use headers from this directory
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
# include_directories($ENV{HOME}/local/include)
# link_directories($ENV{HOME}/local/lib)
# # ============ ZLIB_LIBRARY  ============
# find_library(ZLIB_LIBRARY z)
# # ============ pybind11 (for Python module) ============
# include(FetchContent)
# FetchContent_Declare(
#   pybind11
#   GIT_REPOSITORY https://github.com/pybind/pybind11.git
#   GIT_TAG        v2.12.0
# )
# FetchContent_MakeAvailable(pybind11)

# # Shared source files used by both Python and CLI/test
# set(SharedSources
#     src/ReadFiles.cpp
#     src/ReadBGEN.cpp
# )

# # ============ Python module (pybind11) ============
# add_library(MyGen MODULE
#     src/PyInterface.cpp
#     src/RunPipeline.cpp
#     ${SharedSources}
# )
# target_link_libraries(MyGen PRIVATE
#     pybind11::module
#     ${ZLIB_LIBRARY}
#     zstd
#     deflate
# )
# # Optional: rename output to suppress 'lib' prefix for Python
# set_target_properties(MyGen PROPERTIES PREFIX "")
# # set_target_properties(MyGen PROPERTIES PREFIX "" SUFFIX ".so")
# # ============ CLI/Library for tests ============

# set(Main
#     src/main.cpp
#     src/RunPipeline.cpp
#     ${SharedSources}
# )


# add_executable(runMain ${Main})
# target_link_libraries(runMain
#     PRIVATE zstd deflate ${ZLIB_LIBRARY})

###################################################################################

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


# Set compiler optimization flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")