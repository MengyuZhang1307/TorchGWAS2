cmake_minimum_required(VERSION 3.25)
project(ReadBgen)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use headers from this directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib)

# ============ pybind11 (for Python module) ============
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.12.0
)
FetchContent_MakeAvailable(pybind11)

# Shared source files used by both Python and CLI/test
set(SharedSources
    ReadFiles.cpp
    ReadBGEN.cpp
)

# ============ Python module (pybind11) ============
add_library(MyGen MODULE
    PyInterface.cpp
    RunPipeline.cpp
    ${SharedSources}
)
target_link_libraries(MyGen PRIVATE
    pybind11::module
    zstd
    deflate
)
# Optional: rename output to suppress 'lib' prefix for Python
set_target_properties(MyGen PROPERTIES PREFIX "")
# set_target_properties(MyGen PROPERTIES PREFIX "" SUFFIX ".so")
# ============ CLI/Library for tests ============

add_library(MyCode STATIC
    main.cpp          # You can include this only if you want CLI logic here
    ${SharedSources}
)

# ============ GoogleTest ============

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(runTests
    test/test_bgen.cpp)
target_link_libraries(runTests 
    PRIVATE gtest_main MyCode zstd deflate)
include(GoogleTest)
gtest_discover_tests(runTests)
