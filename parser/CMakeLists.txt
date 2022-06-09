cmake_minimum_required(VERSION 3.0.0)
project(parser VERSION 0.1.0)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(include)
file(GLOB srcs src/*.cpp)
add_executable(parser main.cpp ${srcs})