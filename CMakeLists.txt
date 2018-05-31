cmake_minimum_required(VERSION 3.10)
project(HElib-basic)

set(CMAKE_CXX_STANDARD 11)

# Load in specified order
add_subdirectory(3rdparty)
include_directories(include)
add_subdirectory(source)
add_subdirectory(test)

enable_testing()