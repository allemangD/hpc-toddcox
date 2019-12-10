cmake_minimum_required(VERSION 2.8)
project(hpc)

add_subdirectory(common)

add_subdirectory(cpu-slo)
add_subdirectory(cpu-opt)

add_subdirectory(gpu-slo)
