#pragma once

#include "util.h"

#include <iostream>

/*
 * Order 4*res*res
 */
Coxeter torus(int res) {
    return make_coxeter(4, {
        {0, 1, res},
        {2, 3, res},
    });
}

Coxeter hypercube(int dim) {
    std::vector<Rel> rels;
    rels.push_back({0, 1, 4});
    for (int d = 2; d < dim; d++) {
        rels.push_back({d-1, d, 3});
    }
    return make_coxeter(dim, rels);
}

/*
 * Order 14,400
 */
Coxeter H4() {
    return make_coxeter(4, {
        {0, 1, 5},
        {1, 2, 3},
        {2, 3, 3},
    });
}

/*
 * Order 51,840
 */
Coxeter E6() {
    return make_coxeter(6, {
        {0, 1, 3},
        {1, 2, 3},
        {2, 3, 3},
        {2, 4, 3},
        {4, 5, 3},
    });
}

/*
 * Order 2,903,040
 */
Coxeter E7() {
    return make_coxeter(7, {
        {0, 1, 3},
        {1, 2, 3},
        {2, 3, 3},
        {2, 4, 3},
        {4, 5, 3},
        {5, 6, 3},
    });
}

/*
 * Order 696,729,600
 */
Coxeter E8() {
    return make_coxeter(8, {
        {0, 1, 3},
        {1, 2, 3},
        {2, 3, 3},
        {2, 4, 3},
        {4, 5, 3},
        {5, 6, 3},
        {6, 7, 3},
    });
}

/*
 * returns coxeter group based on the arguments
 * prints out type and arguments, without an endline
 */
Coxeter proc_args(int argc, const char* argv[]) {
    if (argc < 2) {
        std::cerr << "missing type argument." << std::endl;
        exit(EXIT_FAILURE);
    }

    int type = std::strtol(argv[1], nullptr, 10);
    std::cout << type << ',';

    int arg;
    switch (type) {
    case 0:
        if (argc < 3) {
            std::cerr << "Must provide a size for torus!" << std::endl;
            exit(EXIT_FAILURE);
        }
        arg = std::strtol(argv[2], nullptr, 10);
        std::cout << arg << ',';
        return torus(arg);
    case 1:
        std::cout << -1 << ',';
        return H4();
    case 2:
        std::cout << -1 << ',';
        return E6();
    case 3:
        std::cout << -1 << ',';
        return E7();
    case 4:
        std::cout << -1 << ',';
        return E8();
    case 5:
        if (argc < 3) {
            std::cerr << "Must provide a dimension for hypercube!" << std::endl;
            exit(EXIT_FAILURE);
        }
        arg = std::strtol(argv[2], nullptr, 10);
        std::cout << arg << ',';
        return hypercube(arg);
    }

    std::cerr << "Not a valid type!" << std::endl;
    exit(EXIT_FAILURE);
}

