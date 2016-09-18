cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

project (elektron)

# Version
set (ELEKTRON_VERSION_MAJOR 0)
set (ELEKTRON_VERSION_MINOR 1)

add_subdirectory(math)
add_subdirectory(cv)
add_subdirectory(filters)

add_executable(elektron main.cpp)
target_link_libraries (elektron math)
target_link_libraries (elektron cv)
target_link_libraries (elektron filters)
