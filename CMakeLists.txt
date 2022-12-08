cmake_minimum_required(VERSION 3.5)
project(SVClusterer)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.16/include)

add_executable(SVClusterer clusterer.cpp)

find_library(HTS_LIB hts PATHS htslib-1.16/lib NO_DEFAULT_PATH)
target_link_libraries(SVClusterer "${HTS_LIB}")
