#include <cstdlib>

#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>

#include <omp.h>

using Gen=int;
using Gens=std::vector<Gen>;
using Table=std::vector<Gens>;

using Cos=int;

const size_t ALIGN_SIZE=64;
const size_t GENS_PER_LINE=ALIGN_SIZE/sizeof(Gen);

struct Coxeter {
    Gen *gen[2];
    int *size; // multiplicity * 2
    int ngens;
    int nrels;
    Coxeter(int ngens): ngens(ngens){
        nrels = (ngens*(ngens-1))>>1;
        int allocsize = ((nrels-1)/(GENS_PER_LINE) + 1)*ALIGN_SIZE;
        int err = 0;
        err |= posix_memalign((void**)&(gen[0]),ALIGN_SIZE,allocsize);
        err |= posix_memalign((void**)&(gen[1]),ALIGN_SIZE,allocsize);
        err |= posix_memalign((void**)&(size),ALIGN_SIZE,allocsize);
        if (err != 0) {
            std::cerr << "Error allocating Coxeter!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    void clean() {
        free(gen[0]);
        free(gen[1]);
        free(size);
    }
};

struct Mult {
    int from, to, multiplicity;
};

Coxeter make_coxeter(int ngens, const std::vector<Mult> &ms) {
    int mults[ngens][ngens];

    for (int i = 0; i < ngens; i++) {
        for (int j = 0; j < ngens; j++) {
            mults[i][j] = 2;
            mults[j][i] = 2;
        }
    }

    for (const auto &m : ms) {
        mults[m.from][m.to] = m.multiplicity;
        mults[m.to][m.from] = m.multiplicity;
    }

    Coxeter c(ngens);
    int k=0;
    for (int i = 0; i < ngens; i++) {
        for (int j = i + 1; j < ngens; j++) {
            int size = mults[i][j]<<1;
            c.gen[0][k] = i;
            c.gen[1][k] = j;
            c.size[k] = size;
            k++;
        }
    }

    return c;
}

/*
 * Order 4*res*res
 */
Coxeter torus(int res) {
    return make_coxeter(4, {
        {0, 1, res},
        {2, 3, res},
    });
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

void pp(const Gens &g, int w) {
    for (const auto &e : g) {
        std::cerr << std::setw(w) << e << " ";
    }
    std::cerr << std::endl;
}

void pp(const Table &t) {
    std::cerr << "| table:" << std::endl;
    int w = 3;

    for (size_t i = 0; i < t.size(); ++i) {
        std::cerr << std::setw(w) << i << " | ";
        pp(t[i], w);
    }
}

void add_row(const int ngens, const Coxeter &cox,
    Table &cosets, std::vector<Table> &reltables,
    Table &starts, Table &ends) {

    int C = cosets.size();

    cosets.emplace_back(ngens, -1);

    for (unsigned int i = 0; i < cox.nrels; ++i) {
        auto &table = reltables[i];

        unsigned int R = cox.size[i];

        table.emplace_back(R + 1, -1);
        table[C][0] = C;
        table[C][R] = C;

        starts[i].push_back(0);
        ends[i].push_back(R);
    }
}

int add_coset(const int ngens, const Coxeter &cox,
    Table &cosets, std::vector<Table> &reltables,
    Table &starts, Table &ends,
    int coset_scan_hint) {

    int C = cosets.size();

    for (int c = coset_scan_hint; c < C; ++c) {
        std::vector<int> &row = cosets[c];
        for (int g = 0; g < ngens; ++g) {
            if (row[g] == -1) {
                row[g] = C;
                add_row(ngens, cox, cosets, reltables, starts, ends);
                cosets[C][g] = c;
                return c;
            }
        }
    }

    return -1;
}


/**
 * learn until it can't
 */
void learn(Table &coset, const Coxeter &cox,
    std::vector<Table> &reltables, Table &starts, Table &ends) {

    unsigned int nrels = cox.nrels;

    while (true) {
        bool complete = true;

        Gen gens[2];
#pragma omp parallel for schedule(static, 1) reduction(&:complete) private(gens)
        for (unsigned int r = 0; r < nrels; ++r) {
            auto &table = reltables[r];
            gens[0] = cox.gen[0][r];
            gens[1] = cox.gen[1][r];

            for (unsigned int c = 0; c < table.size(); c++) {
                auto &row = table[c];
                auto s = starts[r][c];
                auto e = ends[r][c];

                if (s == e - 1) continue;

                while (row[s + 1] == -1) {
                    const int &lookup = coset[row[s]][gens[s&1]];
                    if (lookup < 0) break;

                    s++;
                    row[s] = lookup;
                }

                while (row[e - 1] == -1) {
                    const int &lookup = coset[row[e]][gens[(e - 1)&1]];
                    if (lookup < 0) break;

                    e--;
                    row[e] = lookup;
                }

                if (s == e - 1) {
                    complete = false;

                    const int &gen = gens[s&1];
                    coset[row[s]][gen] = row[e];
                    coset[row[e]][gen] = row[s];
                }

                starts[r][c] = s;
                ends[r][c] = e;
            }
        }

        if (complete) break;
    }
}

Table solve_tc(int ngens, const Gens &subgens, const Coxeter &cox) {
    Table cosets;
    std::vector<Table> reltables(cox.nrels);

    // storing progress for each relation table row
    Table starts(cox.nrels); // [rel_table][coset]
    Table ends(cox.nrels);

    // set up initial coset
    add_row(ngens, cox, cosets, reltables, starts, ends);
    for (const auto &gen : subgens) {
        cosets[0][gen] = 0;
    }

    int coset_scan_hint = 0;
    while (coset_scan_hint >= 0) {
        learn(cosets, cox, reltables, starts, ends);
        coset_scan_hint = add_coset(ngens, cox, cosets, reltables, starts, ends, coset_scan_hint);
    }

    return cosets;
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
    }

    std::cerr << "Not a valid type!" << std::endl;
    exit(EXIT_FAILURE);
}

int main(int argc, const char *argv[]) {
    Coxeter cox = proc_args(argc, argv);

    auto s = std::chrono::system_clock::now();
    auto cosets = solve_tc(cox.ngens, {}, cox);
    auto e = std::chrono::system_clock::now();

    std::chrono::duration<float> diff = e - s;
    size_t order = cosets.size();

    // type,arg,ngens,time,order
    std::cout << cox.ngens << ',' << diff.count() << ',' << order << std::endl;

    cox.clean();
    return 0;
}
