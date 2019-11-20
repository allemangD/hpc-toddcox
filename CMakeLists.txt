cmake_minimum_required(VERSION 3.14)
project(hpc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(OpenMP)

add_executable(hpc main.cpp)
target_link_libraries(hpc PUBLIC OpenMP::OpenMP_CXX)
