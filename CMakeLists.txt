# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(subnucleondiffraction CXX C)

set(CMAKE_BUILD_TYPE Release)


find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})


add_subdirectory(src)


set_target_properties(wline PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
