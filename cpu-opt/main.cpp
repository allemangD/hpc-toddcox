#include <cstdlib>
#include <immintrin.h>
#include "aligned_allocator.cpp"

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
    std::vector<Cos,aligned_allocator<Cos,ALIGN_SIZE>> start_cosets;
    std::vector<Cos,aligned_allocator<Cos,ALIGN_SIZE>> end_cosets;
    std::vector<Ind,aligned_allocator<Cos,ALIGN_SIZE>> start_inds;
    std::vector<Ind,aligned_allocator<Cos,ALIGN_SIZE>> end_inds;
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

struct CosetTable {
    std::vector<Cos> table;
    int num_cosets;
    int ngens;
    CosetTable (int ngens): ngens(ngens), num_cosets(0) {}
    void add_row() {
        num_cosets++;
        table.resize(table.size()+ngens, -1);
    }
    inline Cos &operator[](int idx) {
        return table[idx];
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

void add_row(const Coxeter &cox,
    CosetTable &cosets, std::vector<RelTable> &reltables) {

    int C = cosets.num_cosets;

    cosets.add_row();

    for (RelTable &rt : reltables)
        rt.add_row(C);
}

int add_coset(const Coxeter &cox,
    CosetTable &cosets, std::vector<RelTable> &reltables,
    int coset_scan_hint) {

    const int C = cosets.num_cosets;
    const int ngens = cox.ngens;

    for (int c = coset_scan_hint; c < C; ++c) {
        int idx = c*ngens;
        for (int g = 0; g < ngens; ++g) {
            if (cosets[idx] == -1) {
                cosets[idx] = C;
                add_row(cox, cosets, reltables);
                cosets[C*ngens+g] = c;
                return c;
            }
            idx++;
        }
    }

    return -1;
}

#ifndef __AVX2__

/**
 * learn until it can't
 */
void learn(const Coxeter &cox, CosetTable &cosets,
    std::vector<RelTable> &reltables) {

    const int nrels = cox.nrels;
    const int ngens = cox.ngens;

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
                    const int lookup = cosets[s_c*ngens + gens[s_i&1]];
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
                    const int lookup = cosets[e_c*ngens + gens[e_i&1]];
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

                    const int gen = gens[s_i&1];
                    cosets[s_c*ngens + gen] = e_c;
                    cosets[e_c*ngens + gen] = s_c;

                    table.rem_row(c);
                    c--;
                }
            }
        }

        if (complete) break;
    }
}

#else

const auto mask = _mm256_set1_epi32(1);
inline void compute_lookups(Ind* inds, Cos* coss, Gen* gens, Cos* base, __m256i step, Cos* tar) {
    _mm256_store_si256(
        (__m256i*)tar,
        _mm256_i32gather_epi32(
            (int*)base,
            _mm256_add_epi32(
                _mm256_mullo_epi32(
                    _mm256_load_si256((__m256i*)coss),
                    step
                ),
                _mm256_i32gather_epi32(
                    gens,
                    _mm256_and_si256(
                        _mm256_load_si256((__m256i*)inds),
                        mask
                    ),
                    4
                )
            ),
            4
        )
    );
}

/**
 * learn until it can't
 */
void learn(const Coxeter &cox, CosetTable &cosets,
    std::vector<RelTable> &reltables) {

    const int nrels = cox.nrels;
    const int ngens = cox.ngens;
    const auto step = _mm256_set1_epi32(ngens);

    while (true) {
        bool complete = true;

        alignas(32) Gen gens[2];
        alignas(32) Cos lookups[8];
        alignas(32) Cos init_cosets[8];
        alignas(32) Ind start_inds[8];
        alignas(32) Ind end_inds[8];
        alignas(32) Cos start_cosets[8];
        alignas(32) Cos end_cosets[8];
#pragma omp parallel for schedule(static, 1) reduction(&:complete) private(gens, lookups, init_cosets, start_inds, end_inds, start_cosets, end_cosets)
        for (unsigned int r = 0; r < nrels; ++r) { 
            auto &table = reltables[r];
            gens[0] = table.gen[0];
            gens[1] = table.gen[1];
            
            unsigned int c;
            bool redo_cval;
            for (c = 0; c < ((table.num_rows>>3)<<3); c+=8) {
                redo_cval = false;
                for (int c_ = 0; c_ < 8; c_++) {
                    init_cosets[c_] = table.init_cosets[c+c_];
                    start_inds[c_] = table.start_inds[c+c_];
                    end_inds[c_] = table.end_inds[c+c_];
                    start_cosets[c_] = table.start_cosets[c+c_];
                    end_cosets[c_] = table.end_cosets[c+c_];
                }

                bool startdone = false;
                bool need_reload = false;
                int reload_idx;
                int idx, lookup;
                while (!startdone) {
                    startdone = true;
                    if (need_reload and idx < table.num_rows) {
                        init_cosets[reload_idx] = table.init_cosets[idx];
                        start_inds[reload_idx] = table.start_inds[idx];
                        end_inds[reload_idx] = table.end_inds[idx];
                        start_cosets[reload_idx] = table.start_cosets[idx];
                        end_cosets[reload_idx] = table.end_cosets[idx];
                    }
                    compute_lookups(start_inds, start_cosets, gens, &(cosets[0]), step, lookups);
                    for (int c_ = 0; c_ < 8; c_++) {
                        lookup = lookups[c_];
                        if (start_inds[c_] < end_inds[c_] and lookup >= 0) {
                            start_inds[c_]++;
                            start_cosets[c_] = lookups[c_];
                            startdone = false;
                            if (lookup > init_cosets[c_]) {
                                int idx = table.coset_poss[lookup];
                                if (idx >= 0) {
                                    table.rem_row(idx);
                                    if ( (idx>>3)<<3 == c ) {
                                        need_reload = true;
                                        reload_idx = (idx & 8);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                bool enddone = false;
                while (!enddone) {
                    enddone = true;
                    if (need_reload and idx < table.num_rows) {
                        init_cosets[reload_idx] = table.init_cosets[idx];
                        start_inds[reload_idx] = table.start_inds[idx];
                        end_inds[reload_idx] = table.end_inds[idx];
                        start_cosets[reload_idx] = table.start_cosets[idx];
                        end_cosets[reload_idx] = table.end_cosets[idx];
                    }
                    compute_lookups(end_inds, end_cosets, gens, &(cosets[0]), step, lookups);
                    for (int c_ = 0; c_ < 8; c_++) {
                        lookup = lookups[c_];
                        if (start_inds[c_] < end_inds[c_] and lookup >= 0) {
                            end_inds[c_]--;
                            end_cosets[c_] = lookups[c_];
                            enddone = false;
                            if (lookup > init_cosets[c_]) {
                                int idx = table.coset_poss[lookup];
                                if (idx >= 0) {
                                    table.rem_row(idx);
                                    if ( (idx>>3)<<3 == c ) {
                                        need_reload = true;
                                        reload_idx = (idx & 8);
                                        redo_cval = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                for (int c_ = 0; c_ < 8 and c+c_ < table.num_rows; c_++) {
                    table.start_inds[c+c_] = start_inds[c_];
                    table.end_inds[c+c_] = end_inds[c_];
                    table.start_cosets[c+c_] = start_cosets[c_];
                    table.end_cosets[c+c_] = end_cosets[c_];
                }

                for (int c_ = 8; c_ >= 0; c_--) {
                    if (c+c_ < table.num_rows)
                        continue;
                    Ind s_i = start_inds[c_];
                    if (start_inds[s_i] == end_inds[c_]) {
                        complete = false;

                        const Gen gen = gens[s_i&1];
                        Cos s_c = start_cosets[c_];
                        Cos e_c = end_cosets[c_];
                        cosets[s_c*ngens + gen] = e_c;
                        cosets[e_c*ngens + gen] = s_c;
                        
                        table.rem_row(c+c_);
                        redo_cval = true;
                    }
                }
                
                if (redo_cval)
                    c-=8;
            }

            for (; c < table.num_rows; c++) {
                auto s_i = table.start_inds[c];
                auto e_i = table.end_inds[c];
                auto s_c = table.start_cosets[c];
                auto e_c = table.end_cosets[c];
                auto i_c = table.init_cosets[c];

                while (s_i < e_i) {
                    const Cos lookup = cosets[s_c*ngens + gens[s_i&1]];
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
                    const Cos lookup = cosets[e_c*ngens + gens[e_i&1]];
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

                    const Gen gen = gens[s_i&1];
                    cosets[s_c*ngens + gen] = e_c;
                    cosets[e_c*ngens + gen] = s_c;

                    table.rem_row(c);
                    c--;
                }
            }
        }

        if (complete) break;
    }
}

#endif

CosetTable solve_tc(const Coxeter &cox, const Gens &subgens) {
    CosetTable cosets(cox.ngens);
    std::vector<RelTable> reltables;

    for (int i=0; i<cox.nrels; i++) {
        reltables.emplace_back(cox.gen[0][i], cox.gen[1][i], cox.size[i]-1);
    }

    // set up initial coset
    add_row(cox, cosets, reltables);
    for (const auto &gen : subgens) {
        cosets[gen] = 0;
    }

    int coset_scan_hint = 0;
    char a;
    while (coset_scan_hint >= 0) {
        learn(cox, cosets, reltables);
        coset_scan_hint = add_coset(cox, cosets, reltables, coset_scan_hint);
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
    auto cosets = solve_tc(cox, {});
    auto e = std::chrono::system_clock::now();

    std::chrono::duration<float> diff = e - s;
    size_t order = cosets.num_cosets;

    // type,arg,ngens,time,order
    std::cout << cox.ngens << ',' << diff.count() << ',' << order << std::endl;

    cox.clean();
    return 0;
}
