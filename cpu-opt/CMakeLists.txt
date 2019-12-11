set(CMAKE_CXX_FLAGS "-std=gnu++11 -fopenmp -O3")

add_executable(cpu-opt main.cpp solver.hpp)
