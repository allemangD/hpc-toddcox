#include <iostream>
#include <vector>

using Gens=std::vector<int>;
using Table=std::vector<Gens>;

void pp(const Gens &g) {
    for (const auto &e : g) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
}

void pp(const Table &t) {
    std::cout << "| table:" << std::endl;
    for (const auto &g : t) {
        std::cout << "| ";
        pp(g);
    }
}

void new_coset(const int ngens, const std::vector<Gens> &rels,
    Table &cosets, std::vector<Table> &reltables,
    Table &starts, Table &ends) {

    int C = cosets.size();

    cosets.emplace_back(ngens, -1);

    for (size_t i = 0; i < rels.size(); ++i) {
        auto &table = reltables[i];

        unsigned int R = rels[i].size();

        table.emplace_back(R + 1, -1);
        table[C][0] = C;
        table[C][R] = C;

        starts[i].push_back(0);
        ends[i].push_back(R);
    }
}

Table solve_tc(int ngens, const Gens &sub, const std::vector<Gens> &rels) {
    Table cosets;
    std::vector<Table> reltables(rels.size());

    // storing progress for each relation table row
    Table starts(rels.size()); // [rel_table][coset]
    Table ends(rels.size());

    new_coset(ngens, rels, cosets, reltables, starts, ends);

    return cosets;
}

int main() {
    auto cosets = solve_tc(3, {0, 1}, {
        {0, 1, 0, 1, 0, 1, 0, 1},
        {1, 2, 1, 2, 1, 2},
        {0, 2, 0, 2},
    });

    pp(cosets);

    return 0;
}
