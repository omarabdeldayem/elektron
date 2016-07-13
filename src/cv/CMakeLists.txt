file(GLOB CV_SOURCES *.cpp)
add_library (cv ${CV_SOURCES})

include_directories (cv ${CMAKE_CURRENT_SOURCE_DIR})
