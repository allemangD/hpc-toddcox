#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>

#include <omp.h>

using Gens=std::vector<int>;
using Table=std::vector<Gens>;

void pp(const Gens &g, int w) {
    for (const auto &e : g) {
        std::cout << std::setw(w) << e << " ";
    }
    std::cout << std::endl;
}

void pp(const Table &t) {
    std::cout << "| table:" << std::endl;
    int w = 3;

    for (size_t i = 0; i < t.size(); ++i) {
        std::cout << std::setw(w) << i << " | ";
        pp(t[i], w);
    }
}

void add_row(const int ngens, const std::vector<Gens> &rels,
    Table &cosets, std::vector<Table> &reltables,
    Table &starts, Table &ends) {

    int C = cosets.size();

    cosets.emplace_back(ngens, -1);

    for (unsigned int i = 0; i < rels.size(); ++i) {
        auto &table = reltables[i];

        unsigned int R = rels[i].size();

        table.emplace_back(R + 1, -1);
        table[C][0] = C;
        table[C][R] = C;

        starts[i].push_back(0);
        ends[i].push_back(R);
    }
}

int add_coset(const int ngens, const std::vector<Gens> &rels,
    Table &cosets, std::vector<Table> &reltables,
    Table &starts, Table &ends,
    int coset_scan_hint) {

    int C = cosets.size();

    for (int c = coset_scan_hint; c < C; ++c) {
        std::vector<int> &row = cosets[c];
        for (int g = 0; g < ngens; ++g) {
            if (row[g] == -1) {
                row[g] = C;
                add_row(ngens, rels, cosets, reltables, starts, ends);
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
void learn(Table &coset, const std::vector<Gens> &rels,
    std::vector<Table> &reltables, Table &starts, Table &ends) {

    unsigned int nrels = rels.size();

    while (true) {
        bool complete = true;

#pragma omp parallel for schedule(static, 1) reduction(&:complete)
        for (unsigned int r = 0; r < nrels; ++r) {
            auto &table = reltables[r];
            const auto &rel = rels[r];

            for (unsigned int c = 0; c < table.size(); c++) {
                auto &row = table[c];
                auto s = starts[r][c];
                auto e = ends[r][c];

                if (s == e - 1) continue;

                while (row[s + 1] == -1) {
                    const int &lookup = coset[row[s]][rel[s]];
                    if (lookup < 0) break;

                    s++;
                    row[s] = lookup;
                }

                while (row[e - 1] == -1) {
                    const int &lookup = coset[row[e]][rel[e - 1]];
                    if (lookup < 0) break;

                    e--;
                    row[e] = lookup;
                }

                if (s == e - 1) {
                    complete = false;

                    const int &gen = rel[s];
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

Table solve_tc(int ngens, const Gens &subgens, const std::vector<Gens> &rels) {
    Table cosets;
    std::vector<Table> reltables(rels.size());

    // storing progress for each relation table row
    Table starts(rels.size()); // [rel_table][coset]
    Table ends(rels.size());

    // set up initial coset
    add_row(ngens, rels, cosets, reltables, starts, ends);
    for (const auto &gen : subgens) {
        cosets[0][gen] = 0;
    }

    int coset_scan_hint = 0;
    while (coset_scan_hint >= 0) {
        learn(cosets, rels, reltables, starts, ends);
        coset_scan_hint = add_coset(ngens, rels, cosets, reltables, starts, ends, coset_scan_hint);
    }

    return cosets;
}

struct Mult {
    int from, to, multiplicity;
};

Table mults(std::vector<Mult> ms) {
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

Table ezmults(int ngens, std::vector<Mult> ms) {
   bool table[ngens][ngens];

   for (int i=0; i<ngens; i++) {
      for (int j=0; j<ngens; j++) {
         table[i][j] = false;
      }
   }

   for (const auto &m : ms) {
      table[m.from][m.to] = true;
      table[m.to][m.from] = true;
   }

   Table res = mults(ms);
   for (int i=0; i<ngens; i++) {
      for (int j=i+1; j<ngens; j++) {
         if (!table[i][j]) {
            res.push_back({i,j,i,j});
         }
      }
   }
   return res;
}

/*
 * Order 4*res*res
 */
std::pair<Table,int> torus(int res) {
    return std::make_pair(ezmults(4,{
        {0,1,res},
        {2,3,res},
    }),4);
}

/*
 * Order 14,400
 */
std::pair<Table,int> H4() {
   return std::make_pair(ezmults(4,{
      {0,1,5},
      {1,2,3},
      {2,3,3},
   }),4);
}

/*
 * Order 51,840
 */
std::pair<Table,int> E6() {
   return std::make_pair(ezmults(6,{
      {0,1,3},
      {1,2,3},
      {2,3,3},
      {2,4,3},
      {4,5,3},
   }),6);
}

/*
 * Order 2,903,040
 */
std::pair<Table,int> E7() {
   return std::make_pair(ezmults(7,{
      {0,1,3},
      {1,2,3},
      {2,3,3},
      {2,4,3},
      {4,5,3},
      {5,6,3},
   }),7);
}

/*
 * Order 696,729,600
 */
std::pair<Table,int> E8() {
   return std::make_pair(ezmults(8,{
      {0,1,3},
      {1,2,3},
      {2,3,3},
      {2,4,3},
      {4,5,3},
      {5,6,3},
      {6,7,3},
   }),8);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "REEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl;
        std::cerr << "gimme more ~~tendies~~ arguments" << std::endl;
        return 1;
    }

    int type = std::strtol(argv[1], nullptr, 10);
    int arg = 0;
    std::pair<Table,int> res;
    switch (type) {
      case 0:
         if (argc < 3) {
            std::cerr << "Must provide a size for torus!" << std::endl;
            return 2;
         }
         arg = std::strtol(argv[2], nullptr, 10);
         res = torus(arg);
         break;
      case 1:
         res = H4();
         break;
      case 2:
         res = E6();
         break;
      case 3:
         res = E7();
         break;
      case 4:
         res = E8();
         break;
      default:
         std::cerr << "Not a valid type!" << std::endl;
         return 3;
    }

    const Table &rels = res.first;
    const int ngens = res.second;

    auto s = std::chrono::system_clock::now();
    auto cosets = solve_tc(ngens, {}, rels);
    auto e = std::chrono::system_clock::now();

    std::chrono::duration<float> diff = e - s;
    size_t order = cosets.size();

    // ngens,type,arg,time,order
    std::cout << ngens << ',' << type << ',' << arg  << ',' << diff.count() << ',' << order << std::endl;

    return 0;
}
