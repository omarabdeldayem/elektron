file(GLOB MATH_SOURCES *.cpp)
list(REMOVE_ITEM MATH_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
add_library(math ${MATH_SOURCES})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
