file(GLOB STATS_SOURCES *.cpp)
add_library (stats ${STATS_SOURCES})

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
