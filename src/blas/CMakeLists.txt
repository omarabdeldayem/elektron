file(GLOB BLAS_SOURCES *.cpp)
add_library (blas ${BLAS_SOURCES})

target_include_directories (${CMAKE_CURRENT_SOURCE_DIR})
