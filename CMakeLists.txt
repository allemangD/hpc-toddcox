cmake_minimum_required(VERSION 2.8)
project(hpc)

set(CMAKE_CXX_FLAGS "-std=gnu++11 -fopenmp -O3")

add_executable(hpc main.cpp)
