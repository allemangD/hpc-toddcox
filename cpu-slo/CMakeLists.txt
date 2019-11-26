set(CMAKE_CXX_FLAGS "-std=gnu++11 -fopenmp -O3")

add_executable(cpu-slo main.cpp)
target_link_libraries(cpu-slo PRIVATE common)
