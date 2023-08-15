cmake_minimum_required(VERSION 2.8)
project(SVClusterer)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.16/include)

add_executable(clusterer clusterer.cpp)

add_executable(compare-ins libs/ssw.c libs/ssw_cpp.cpp compare-ins.cpp)
add_executable(compare-del libs/ssw.c libs/ssw_cpp.cpp compare-del.cpp)

#add_executable(annotate_svs annotate_svs.cpp)

find_library(HTS_LIB hts PATHS htslib-1.16/lib NO_DEFAULT_PATH)
target_link_libraries(clusterer "${HTS_LIB}")
target_link_libraries(compare-del "${HTS_LIB}")
target_link_libraries(compare-ins "${HTS_LIB}")
#target_link_libraries(annotate_svs "${HTS_LIB}")
