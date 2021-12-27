cmake_minimum_required(VERSION 3.21)
project(algorithms_template)

set(CMAKE_CXX_STANDARD 17)

add_executable(algorithms_template LIS.cpp choose.cpp TSP.cpp articulation_and_bridges.cpp kosaraju.cpp prim.cpp prim2.cpp modifieddijkstras.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -D_GLIBCXX_DEBUG -fsanitize=address -fsanitize=undefined")
