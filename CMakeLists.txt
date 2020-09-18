cmake_minimum_required (VERSION 3.4)
project (kDeduper C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -lgomp -lrt -fopenmp -W -Wall -pedantic -O3 -Ofast -march=native")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")

# Other libs
include_directories(libs/kseq/)
include_directories(libs/parallel_hashmap/)
include_directories(includes)

# Executables
add_executable (kDeduper kDeduper.cpp src/kmer.cpp)
target_link_libraries (kDeduper z)