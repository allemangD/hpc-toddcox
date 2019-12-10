#pragma once

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <iostream>

template<class T>
std::ostream &operator<<(std::ostream &o, const thrust::host_vector<T> &vec) {
    if (vec.size() == 0 || vec.size() > 15) 
        return o << "host_vector{size=" << vec.size() << "}";

    o << "[";
    
    for (int i = 0; i < vec.size() - 1; i++) o << vec[i] << ", ";
    
    if (vec.size() > 0) o << vec[vec.size() - 1];
    
    o << "]";

    return o;
}

template<class T>
std::ostream &operator<<(std::ostream &o, const thrust::device_vector<T> &vec) {
    return o << "device_vector{size=" << vec.size() << "}";
}

