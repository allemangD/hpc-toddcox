cmake_minimum_required(VERSION 2.8)
project(hpc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)

add_executable(hpc main.cpp)
target_link_libraries(hpc PUBLIC OpenMP::OpenMP_CXX)
