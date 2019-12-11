cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(main LANGUAGES CXX CUDA)

add_executable(main main.cu)
target_compile_features(main PUBLIC cxx_std_11)
set_target_properties(main PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
