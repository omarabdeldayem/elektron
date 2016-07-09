cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

project (NEUROLIB)

# Version
set (NEUROLIB_VERSION_MAJOR 0)
set (NEUROLIB_VERSION_MINOR 1)

add_subdirectory(blas)
#add_subdirectory(cv)

add_executable(neurolib Main.cpp)
target_link_libraries (neurolib blas)
