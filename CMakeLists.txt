cmake_minimum_required(VERSION 3.25)
project(ReadBgen)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(mycode main.cpp
            ReadFiles.cpp
            ReadBGEN.cpp)

include_directories($ENV{HOME}/local/include)
link_directories($ENV{HOME}/local/lib)

# Download and link GoogleTest############################################################
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()


# Create test executable
add_executable(runTests test/test_bgen.cpp)
target_link_libraries(runTests gtest_main mycode zstd z deflate)

include(GoogleTest)
gtest_discover_tests(runTests)