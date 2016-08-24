file(GLOB BLAS_SOURCES *.cpp)
list(REMOVE_ITEM BLAS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Matrix.hpp)
add_library(blas ${BLAS_SOURCES})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
