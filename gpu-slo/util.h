#pragma once

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>
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

struct Rel {
    int gens[2];
    int mul;
};

struct Coxeter {
    int ngens;
    std::vector<Rel> rels;
};

Coxeter make_coxeter(int ngens, const std::vector<Rel> &rels) {
    int mults[ngens][ngens];

    for (int i = 0; i < ngens; i++) {
        for (int j = 0; j < ngens; j++) {
            mults[i][j] = 2;
            mults[j][i] = 2;
        }
    }

    for (const auto &r : rels) {
        mults[r.gens[0]][r.gens[1]] = r.mul;
        mults[r.gens[1]][r.gens[0]] = r.mul;
    }

    std::vector<Rel> res;
    
    for (int i = 0; i < ngens; i++) {
        for (int j = i + 1; j < ngens; j++) {
            res.push_back({i, j, mults[i][j]});
        }
    }

    return {ngens, res};
}
