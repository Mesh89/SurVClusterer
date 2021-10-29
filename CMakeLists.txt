cmake_minimum_required(VERSION 3.15)
project(SVClusterer)

set(CMAKE_CXX_STANDARD 14)

include_directories(htslib-1.13/include)

add_executable(SVClusterer clusterer.cpp)

find_library(HTS_LIB hts PATHS htslib-1.13/lib NO_DEFAULT_PATH)
target_link_libraries(SVClusterer "${HTS_LIB}")
