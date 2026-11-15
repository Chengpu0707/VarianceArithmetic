/*
Micro-benchmark for Matrix::calc() at sizes 6, 8, 10. Not part of the test
suite — build/run on demand:
    g++ -std=c++23 -O2 -pthread -o BenchMatrix.exe BenchMatrix.cpp -lstdc++exp
    .\BenchMatrix.exe
*/
#if __cplusplus >= 201103L
#include <chrono>
#else
#include <ctime>
#endif
#include <iostream>
#include "Matrix.h"

using namespace var_dbl;

// Elapsed milliseconds since the previous call (wall clock; CPU clock before C++11).
static double elapsedMs() {
#if __cplusplus >= 201103L
    typedef std::chrono::high_resolution_clock clk;
    static clk::time_point t0 = clk::now();
    const clk::time_point t1 = clk::now();
    const double ms = std::chrono::duration<double>(t1 - t0).count() * 1000.0;
#else
    static std::clock_t t0 = std::clock();
    const std::clock_t t1 = std::clock();
    const double ms = 1000.0 * static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
#endif
    t0 = t1;
    return ms;
}

int main() {
    const size_t sSize[] = {6, 8, 10};
    for (size_t i = 0; i < sizeof(sSize) / sizeof(sSize[0]); ++i) {
        const size_t sz = sSize[i];
        for (int w = 0; w < 3; ++w) (void) Matrix::randomMatrix(sz, 1e-6).determ();
        const int reps = (sz == 10) ? 5 : (sz == 8) ? 60 : 400;
        elapsedMs();
        for (int r = 0; r < reps; ++r) {
            Matrix m = Matrix::randomMatrix(sz, 1e-6);
            VarDbl d = m.determ();
            if (d.value() != d.value()) return 1;
        }
        const double ms = elapsedMs();
        std::cout << "size=" << sz << " reps=" << reps
                  << " total=" << ms << " ms (avg " << (ms / reps) << " ms/matrix)\n";
    }
    return 0;
}
