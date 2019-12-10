#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/for_each.h>

#include <vector>
#include <iostream>

#include "util.h"

// struct Cosets {
//     int width;
//     thrust::device_vector<int> data{};
//     
//     __host__
//     Cosets(int ngens) : width(ngens) {
//     }
//     
//     void add_row() {
//         data.resize(data.size() + width, -1);
//     }
//     
//     thrust::host_vector<int> get_data() {
//         return data;
//     }
// };

struct Row {
    int rel;

    int l, r;

    int from, to;

    __host__ __device__
    Row() : rel(0), l(0), r(0), from(0), to(0) {}

    __device__
    Row(int rel, int cos, int size) {
        l = 0;
        r = size - 1;
        
        from = to = cos;
        
        this->rel = rel;
    }
};

std::ostream &operator<<(std::ostream &o, const Row &r) {
    return o << "Row[" << r.rel << "]{" << r.l << ":" << r.from << "-" << r.to << ":" << r.r << "}";
}

struct Rel {
    int gens[2];
    int mul;
};

struct Solver {
    int ngens;
    int *cosets;
    Rel *rels;
    
    Solver(int ngens,
           thrust::device_vector<int> &cosets,
           thrust::device_vector<Rel> &rels)
        : ngens(ngens),
           cosets(thrust::raw_pointer_cast(cosets.data())),
           rels(thrust::raw_pointer_cast(rels.data())) {
        }
    
    __device__
    void operator()(Row &r) {
        if (r.l + 1 >= r.r) return;
        
        while ((r.r - r.l) > 0) {
            int gen = rels[r.rel].gens[r.l & 1];
            int next = cosets[r.from * ngens + gen];
            if (next < 0) break;
            r.l++;
            r.from = next;
        }

        while ((r.r - r.l) > 0) {
            int gen = rels[r.rel].gens[r.r & 1];
            int next = cosets[r.to * ngens + gen];
            if (next < 0) break;
            r.r--;
            r.to = next;
        }
            
        if (r.r - r.l == 0) {
            int gen = rels[r.rel].gens[r.l & 1];
            cosets[r.from * ngens + gen] = r.to;
            cosets[r.to * ngens + gen] = r.from;
        }
    }
};

struct CosetInitializer {
    int *cosets;

    CosetInitializer(thrust::device_vector<int> &cosets)
        : cosets(thrust::raw_pointer_cast(cosets.data())) {
    }

    __device__
    void operator()(int gen) {
        cosets[gen] = 0;
    }
};

struct RowGen {
    Rel *rels;

    int coset;

    RowGen(int coset, thrust::device_vector<Rel> &rels) 
        : coset(coset),
          rels(thrust::raw_pointer_cast(rels.data())) {}

    __device__
    Row operator()(int rel) {
        return Row(rel, coset, rels[rel].mul * 2);
    }
};

thrust::device_vector<int> solve(
        int ngens,
        thrust::device_vector<int> subs,
        thrust::device_vector<Rel> rels) {
    
    thrust::device_vector<int> cosets;
    cosets.resize(cosets.size() + ngens, -1);
    int lastCoset = 0;

    thrust::for_each(subs.begin(), subs.end(), 
        CosetInitializer(cosets));

    thrust::device_vector<Row> rows;

    thrust::counting_iterator<int> counter(0);

    thrust::device_vector<Row> new_rows(rels.size());
    thrust::transform(counter, counter + rels.size(), new_rows.begin(),
        RowGen(lastCoset, rels));
    rows.insert(rows.begin(), new_rows.begin(), new_rows.end());

    std::cout << rows << std::endl;

    Solver solv(ngens, cosets, rels);

    std::cout << thrust::host_vector<Row>(rows) << std::endl;
    thrust::for_each(rows.begin(), rows.end(), solv);
    std::cout << thrust::host_vector<Row>(rows) << std::endl;

    return cosets;
}


int main(int argc, char* argv[]) {
    int ngens = 4;
    std::vector<Rel> rels = {
        {0, 1, 4},
        {1, 2, 3},
        {2, 3, 3},

        {0, 2, 2},
        {1, 2, 2},
        {1, 3, 2},
    };
    std::vector<int> subs = {1, 3};

    thrust::host_vector<int> cosets = solve(ngens, subs, rels);

    std::cout << cosets << std::endl;

    return 0;
}

