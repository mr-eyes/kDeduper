cmake_minimum_required (VERSION 3.4)
project (kDeduper C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -lgomp -lrt -fopenmp -W -Wall -pedantic -O3 -Ofast -march=native")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")

# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------
add_subdirectory(libs/kProcessor)
set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/libs/kProcessor/include")
include_directories("${PROJECT_SOURCE_DIR}/libs/kProcessor/ThirdParty/MQF/include")
include_directories (${kProcessor_INCLUDE_PATH})
# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

# Other libs
include_directories(libs/kseq/)
include_directories(includes)

# Executables
add_executable (kDeduper kDeduper.cpp)
target_link_libraries (kDeduper kProcessor z)
target_include_directories(kDeduper INTERFACE ${kProcessor_INCLUDE_PATH})

