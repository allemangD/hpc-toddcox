//#include <cstdio>
//#include <cstdlib>
//
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
//#include <thrust/sequence.h>
//
//#define N 50
//
//__global__
//void vector_add(float* out, float* a, float* b, int n) {
//    for(int i = 0; i < n; i++){
//        out[i] = a[i] + b[i];
//    }
//}
//
//int main(){
//    thrust::host_vector<float> a(N);
//    thrust::sequence(a.begin(), a.end());
//
//    thrust::host_vector<float> b(N);
//    thrust::sequence(b.begin(), b.end());
//    thrust::reverse(b.begin(), b.end());
//
//    for (int i = 0; i < N; ++i) {
//        printf("%.1f ", a[i]);
//    } printf("\n");
//
//    for (int i = 0; i < N; ++i) {
//        printf("%.1f ", b[i]);
//    } printf("\n");
//
//    thrust::device_vector<float> aD = a;
//    thrust::device_vector<float> bD = b;
//    thrust::device_vector<float> outD(N);
//
//    vector_add<<<1, 1>>>(
//        thrust::raw_pointer_cast(&outD[0]),
//        thrust::raw_pointer_cast(&aD[0]),
//        thrust::raw_pointer_cast(&bD[0]),
//        N);
//
//    thrust::host_vector<float> out = outD;
//
//    for (int i = 0; i < N; ++i) {
//        printf("%.1f ", out[i]);
//    } printf("\n");
//
//    return 0;
//}

#include <cstdio>
#include <cstdlib>
#include <chrono>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>

void add_proc(int *c, int *a, int *b) {
    *c = *a + *b;
}

void test_proc(){
    int a = 0;
    int b = 1;

    auto start = std::chrono::system_clock::now();

    for (int i = 0; i < 1000000; ++i) {
        add_proc(&a, &a, &b);
    }

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<float, std::micro> diff = end - start;

    printf("proc: %d: 1B in %.3f micro\n", a, diff.count());
}

__global__
void add_gpu(int *c, int *a, int *b) {
    *c = *a + *b;
}

void test_gpu(){
    thrust::device_vector<int> vals(2, 0);
    vals[0] = 0;
    vals[1] = 1;
    printf(" gpu: %d: 1B in %.3f micro\n", vals[0], 0.0f);

    int *a = thrust::raw_pointer_cast(&vals[0]);
    int *b = thrust::raw_pointer_cast(&vals[1]);

    add_gpu<<<1, 1>>>(a, a, a);

    printf(" gpu: %d: 1B in %.3f micro\n", vals[0], 0.0f);
}

int main(int argc, char *argv[]) {
    test_proc();
    test_gpu();
}

