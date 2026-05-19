/*
Micro-benchmark for Matrix::calc() at sizes 6, 8, 10. Not part of the test
suite — build/run on demand:
    g++ -std=c++23 -O2 -pthread -o BenchMatrix.exe BenchMatrix.cpp -lstdc++exp
    .\BenchMatrix.exe
*/
#include <chrono>
#include <iostream>
#include "Matrix.h"

using namespace var_dbl;
using clk = std::chrono::high_resolution_clock;

int main() {
    for (size_t sz : {size_t(6), size_t(8), size_t(10)}) {
        for (int w = 0; w < 3; ++w) (void) Matrix::randomMatrix(sz, 1e-6).determ();
        const int reps = (sz == 10) ? 5 : (sz == 8) ? 60 : 400;
        const auto t0 = clk::now();
        for (int r = 0; r < reps; ++r) {
            Matrix m = Matrix::randomMatrix(sz, 1e-6);
            VarDbl d = m.determ();
            if (d.value() != d.value()) return 1;
        }
        const double ms = std::chrono::duration<double>(clk::now() - t0).count() * 1000.0;
        std::cout << "size=" << sz << " reps=" << reps
                  << " total=" << ms << " ms (avg " << (ms / reps) << " ms/matrix)\n";
    }
    return 0;
}
