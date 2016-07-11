file(GLOB CV_SOURCES *.cpp)
add_library (cv ${CV_SOURCES})

target_include_directories (cv PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
