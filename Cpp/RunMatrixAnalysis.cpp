/*
Driver for Matrix::runMatrixAnalysis — the C++ port of Python/runMatrixAnalysis.py.
Build:
    g++ -std=c++23 -O2 -pthread -o RunMatrixAnalysis.exe RunMatrixAnalysis.cpp -lstdc++exp
Run (defaults: sizes 4..9, 20-noise sweep, 1000 cells, outDir = "Output"):
    .\RunMatrixAnalysis.exe
With explicit args:
    .\RunMatrixAnalysis.exe <minSize> <maxSize> <targetCells> <outDir>
*/
#include <cstdlib>
#include <iostream>
#include <string>
#include "Matrix.h"

using namespace var_dbl;

int main(int argc, char** argv) {
    const size_t minSize     = (argc > 1) ? std::atoi(argv[1]) : 4;
    const size_t maxSize     = (argc > 2) ? std::atoi(argv[2]) : 10;
    const size_t targetCells = (argc > 3) ? std::atoi(argv[3]) : 1000;
    const std::string outDir = (argc > 4) ? argv[4] : "Output";

    std::cout << "runMatrixAnalysis: sizes=[" << minSize << ", " << maxSize
              << ") target=" << targetCells << " outDir=" << outDir << '\n';
    Matrix::runMatrixAnalysis(minSize, maxSize, targetCells,
                              Matrix::defaultNoises(), outDir);
    std::cout << "done\n";
    return 0;
}
