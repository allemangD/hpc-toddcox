#pragma once

#include <types.hpp>

struct Mult {
    int from, to, multiplicity;
};

Table mults(const std::vector<Mult>& ms) {
    Table res;
    for (const auto &m : ms) {
        int N = res.size();
        res.emplace_back(m.multiplicity * 2, m.to);
        for (int i = 0; i < m.multiplicity * 2; i += 2) {
            res[N][i] = m.from;
        }
    }
    return res;
}

std::vector<Mult> ezmults(int ngens, const std::vector<Mult> &ms) {
    bool table[ngens][ngens];

    for (int i = 0; i < ngens; i++) {
        for (int j = 0; j < ngens; j++) {
            table[i][j] = false;
        }
    }

    for (const auto &m : ms) {
        table[m.from][m.to] = true;
        table[m.to][m.from] = true;
    }

    std::vector<Mult> res(ms);

    for (int i = 0; i < ngens; i++) {
        for (int j = i + 1; j < ngens; j++) {
            if (!table[i][j]) {
                res.push_back({i, j, 2});
            }
        }
    }

    return res;
}

/*
 * Order 4*res*res
 */
std::pair<Table, int> torus(int res) {

    return std::make_pair(mults(ezmults(4, {
        {0, 1, res},
        {2, 3, res},
    })), 4);
}

/*
 * Order 14,400
 */
std::pair<Table, int> H4() {
    return std::make_pair(mults(ezmults(4, {
        {0, 1, 5},
        {1, 2, 3},
        {2, 3, 3},
    })), 4);
}

/*
 * Order 51,840
 */
std::pair<Table, int> E6() {
    return std::make_pair(mults(ezmults(6, {
        {0, 1, 3},
        {1, 2, 3},
        {2, 3, 3},
        {2, 4, 3},
        {4, 5, 3},
    })), 6);
}

/*
 * Order 2,903,040
 */
std::pair<Table, int> E7() {
    return std::make_pair(mults(ezmults(7, {
        {0, 1, 3},
        {1, 2, 3},
        {2, 3, 3},
        {2, 4, 3},
        {4, 5, 3},
        {5, 6, 3},
    })), 7);
}

/*
 * Order 696,729,600
 */
std::pair<Table, int> E8() {
    return std::make_pair(mults(ezmults(8, {
        {0, 1, 3},
        {1, 2, 3},
        {2, 3, 3},
        {2, 4, 3},
        {4, 5, 3},
        {5, 6, 3},
        {6, 7, 3},
    })), 8);
}
