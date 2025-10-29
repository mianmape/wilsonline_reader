cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(subnucleondiffraction CXX C)

set(CMAKE_BUILD_TYPE Release)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Only add the src directory
add_subdirectory(src)

