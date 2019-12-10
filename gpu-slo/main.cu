#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/logical.h>

#include <vector>
#include <iostream>

#include "util.h"

struct Row {
    int rel;

    int l, r;

    int from, to;

    bool learning;

    __host__ __device__
    Row() : rel(0), l(0), r(0), from(0), to(0), learning(true) {}

    __device__
    Row(int rel, int cos, int size) {
        l = 0;
        r = size - 1;
        
        from = to = cos;
        
        this->rel = rel;

        learning = true;
    }
};

std::ostream &operator<<(std::ostream &o, const Row &r) {
    return o << "Row[" << r.rel << "]{" << r.l << ":" << r.from << "-" << r.to << ":" << r.r << "}(" << r.learning << ")";
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
        if (r.r - r.l <= 1) {
            r.learning = false;
            return;
        }
        
        while (r.r - r.l > 1) {
            int gen = rels[r.rel].gens[r.l & 1];
            int next = cosets[r.from * ngens + gen];
            if (next < 0) break;
            r.l++;
            r.from = next;
        }

        while (r.r - r.l > 1) {
            int gen = rels[r.rel].gens[r.r & 1];
            int next = cosets[r.to * ngens + gen];
            if (next < 0) break;
            r.r--;
            r.to = next;
        }
            
        if (r.r - r.l <= 1) { 
            int gen = rels[r.rel].gens[r.l & 1];
            cosets[r.from * ngens + gen] = r.to;
            cosets[r.to * ngens + gen] = r.from;

            r.learning = true;
            return;
        }

        r.learning = false;
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

struct RowIncomplete {
    __device__
    bool operator()(Row r) {
        return r.r - r.l > 1;
    }
};

struct Relearn {
    __device__
    void operator()(Row &r) {
        r.learning = true;
    }
};

struct Learning {
    __device__
    bool operator()(Row r) {
        return r.learning;
    }
};

void add_row(
        int ngens,
        thrust::device_vector<int> &cosets) {
    cosets.resize(cosets.size() + ngens, -1);
};

// todo: this part is _real_ slow.
bool add_coset(
        int ngens,
        int *coset,
        int *hint,
        thrust::device_vector<int> &cosets) {
    *coset = cosets.size() / ngens;

    // todo: this part especially.
    while (cosets[*hint] >= 0) {
        *hint = *hint + 1;
        if (*hint >= cosets.size()) 
            return true;
    }
    int from = *hint / ngens;
    int gen = *hint % ngens;
    
    add_row(ngens, cosets);
    
    cosets[*hint] = *coset;
    cosets[*coset * ngens + gen] = from;

    return false;
}

void gen_rows(
        int coset,
        thrust::device_vector<Rel> &rels,
        thrust::device_vector<Row> &rows) {
    rows.resize(rows.size() + rels.size());

    thrust::counting_iterator<int> counter(0);
    thrust::transform(
            thrust::device,
            counter, counter + rels.size(),
            rows.end() - rels.size(),
            RowGen(coset, rels));
}

thrust::device_vector<int> solve(
        int ngens,
        thrust::device_vector<int> subs,
        thrust::device_vector<Rel> rels) {
    
    thrust::device_vector<int> cosets;
    thrust::device_vector<Row> rows;

    add_row(ngens, cosets);
    thrust::for_each(
            thrust::device,
            subs.begin(), subs.end(), 
            CosetInitializer(cosets));

    gen_rows(0, rels, rows);

    int coset = 0;
    int hint = 0;

    while (true) {
        thrust::for_each(
                thrust::device, 
                rows.begin(), 
                rows.end(),
                Relearn());

        Solver solve(ngens, cosets, rels);
        while (true) {
            thrust::for_each(
                    thrust::device,
                    rows.begin(), rows.end(),
                    solve);
            bool r = thrust::any_of(
                    thrust::device,
                    rows.begin(), rows.end(),
                    Learning());
            if (!r) break;
        }

        bool done = add_coset(
                ngens,
                &coset, &hint,
                cosets);

        if (done) break;

        gen_rows(coset, rels, rows);

        auto cut = thrust::partition(
                thrust::device, 
                rows.begin(), rows.end(), 
                RowIncomplete());

        rows.erase(cut, rows.end());
    }

    return cosets;
}


int main(int argc, char* argv[]) {
    // int ngens = 4;
    // std::vector<Rel> rels = {
    //     {0, 1, 3},
    //     {1, 2, 3},
    //     {2, 3, 3},

    //     {0, 2, 2},
    //     {1, 2, 2},
    //     {1, 3, 2},
    // };
    // std::vector<int> subs = {};

    int ngens = 4;
    std::vector<Rel> rels = {
        {0, 1, 4},
        {1, 2, 3},
        {2, 3, 3},

        {0, 2, 2},
        {1, 3, 2},
        {0, 3, 2},
    };
    std::vector<int> subs = {};

    thrust::host_vector<int> cosets = solve(ngens, subs, rels);

    std::cout << cosets.size() / ngens << " cosets" << std::endl;
    for (int c = 0; c < cosets.size(); c += ngens) {
        for (int g = c; g < c + ngens; g++ ) {
            std::cout << cosets[g] << " ";
        } std::cout << std::endl;
    }

    return 0;
}

