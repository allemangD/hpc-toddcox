#include <cstdlib>
#include <immintrin.h>

#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>

#include <omp.h>

using Ind=int;

using Gen=int;
using Gens=std::vector<Gen>;
using Table=std::vector<Gens>;

using Cos=int;

/*
 * COXETER GROUP DEFINITIONS
 */

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
    Mult() {}
    Mult(int from, int to, int multiplicity): from(from), to(to), multiplicity(multiplicity) {}
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

Coxeter hypercube(int dim) {
    std::vector<Mult> mults;
    mults.emplace_back(0,1,4);
    for (int i = 2; i < dim; i++) {
        mults.emplace_back(i-1, i, 3);
    }
    return make_coxeter(dim, mults);
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
 * LEARNING / RelTable DIFINITIONS
 */

struct RelTable {
    std::vector<int> coset_poss;
    std::vector<Cos> init_cosets;
    std::vector<Cos> start_cosets;
    std::vector<Cos> end_cosets;
    std::vector<Ind> start_inds;
    std::vector<Ind> end_inds;
    int num_rows;
    Gen gen[2];
    Ind end_ind;
    RelTable(Gen gen0, Gen gen1, Ind end_ind): end_ind(end_ind), num_rows(0) {
        gen[0] = gen0;
        gen[1] = gen1;
    }
    void add_row(Cos new_coset) {
        coset_poss.push_back(num_rows);
        init_cosets.push_back(new_coset);
        start_cosets.push_back(new_coset);
        end_cosets.push_back(new_coset);
        start_inds.push_back(0);
        end_inds.push_back(end_ind);
        num_rows++;
    }
    void rem_row(int idx) {
        num_rows--;

        coset_poss[init_cosets[num_rows]] = idx;
        coset_poss[init_cosets[idx]] = -1;

        init_cosets[idx] = init_cosets[num_rows];
        init_cosets.pop_back();

        start_cosets[idx] = start_cosets[num_rows];
        start_cosets.pop_back();

        end_cosets[idx] = end_cosets[num_rows];
        end_cosets.pop_back();

        start_inds[idx] = start_inds[num_rows];
        start_inds.pop_back();

        end_inds[idx] = end_inds[num_rows];
        end_inds.pop_back();
    }
};

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
    Table &cosets, std::vector<RelTable> &reltables) {

    int C = cosets.size();

    cosets.emplace_back(ngens, -1);

    for (RelTable &rt : reltables)
        rt.add_row(C);
}

int add_coset(const int ngens, const Coxeter &cox,
    Table &cosets, std::vector<RelTable> &reltables,
    int coset_scan_hint) {

    int C = cosets.size();

    for (int c = coset_scan_hint; c < C; ++c) {
        std::vector<int> &row = cosets[c];
        for (int g = 0; g < ngens; ++g) {
            if (row[g] == -1) {
                row[g] = C;
                add_row(ngens, cox, cosets, reltables);
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
    std::vector<RelTable> &reltables) {

    unsigned int nrels = cox.nrels;

    while (true) {
        bool complete = true;

        Gen gens[2];
#pragma omp parallel for schedule(static, 1) reduction(&:complete) private(gens)
        for (unsigned int r = 0; r < nrels; ++r) { 
            auto &table = reltables[r];
            gens[0] = table.gen[0];
            gens[1] = table.gen[1];

            for (unsigned int c = 0; c < table.num_rows; c++) {
                auto s_i = table.start_inds[c];
                auto e_i = table.end_inds[c];
                auto s_c = table.start_cosets[c];
                auto e_c = table.end_cosets[c];
                auto i_c = table.init_cosets[c];

                while (s_i < e_i) {
                    const int &lookup = coset[s_c][gens[s_i&1]];
                    if (lookup < 0) break;

                    s_i++;
                    s_c = lookup;

                    if (s_c > i_c) {
                        int idx = table.coset_poss[s_c];
                        if (idx >= 0)
                            table.rem_row(idx); 
                    }
                }

                table.start_inds[c] = s_i;
                table.start_cosets[c] = s_c;

                while (s_i < e_i) {
                    const int &lookup = coset[e_c][gens[e_i&1]];
                    if (lookup < 0) break;

                    e_i--;
                    e_c = lookup;

                    if (e_c > i_c) {
                        int idx = table.coset_poss[e_c];
                        if (idx >= 0)
                            table.rem_row(idx);
                    }
                }

                table.end_inds[c] = e_i;
                table.end_cosets[c] = e_c;

                if (s_i == e_i) {
                    complete = false;

                    const int &gen = gens[s_i&1];
                    coset[s_c][gen] = e_c;
                    coset[e_c][gen] = s_c;

                    table.rem_row(c);
                    c--;
                }
            }
        }

        if (complete) break;
    }
}

Table solve_tc(int ngens, const Gens &subgens, const Coxeter &cox) {
    Table cosets;
    std::vector<RelTable> reltables;

    for (int i=0; i<cox.nrels; i++) {
        reltables.emplace_back(cox.gen[0][i], cox.gen[1][i], cox.size[i]-1);
    }

    // set up initial coset
    add_row(ngens, cox, cosets, reltables);
    for (const auto &gen : subgens) {
        cosets[0][gen] = 0;
    }

    int coset_scan_hint = 0;
    char a;
    while (coset_scan_hint >= 0) {
        learn(cosets, cox, reltables);
        coset_scan_hint = add_coset(ngens, cox, cosets, reltables, coset_scan_hint);
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
