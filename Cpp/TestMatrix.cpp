/*
Unit tests for Matrix in Matrix.h: construction, get/set, bottom-up Laplace
determinant (Formula 6.1), adjugate (Formula 6.9 variance), the M·adj(M) = |M|·I
identity, variance propagation, and the random-matrix factory.
*/
#include "Matrix.h"
#include "Test.h"

using namespace var_dbl;


// Write a precise scalar at (r, c).
static void put(Matrix& m, size_t r, size_t c, double val) {
    m.set(r, c, VarDbl(val));
}


void testConstruction() {
    Matrix m(3);
    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 3; ++c) {
            test::assertEqual(m.get(r, c).value(), 0.0);
            test::assertEqual(m.get(r, c).uncertainty(), 0.0);
        }
    }
}

void testGetSet() {
    Matrix m(2);
    put(m, 0, 0, 1.0);  put(m, 0, 1, 2.0);
    put(m, 1, 0, 3.0);  put(m, 1, 1, 4.0);
    test::assertEqual(m.get(0, 0).value(), 1.0);
    test::assertEqual(m.get(0, 1).value(), 2.0);
    test::assertEqual(m.get(1, 0).value(), 3.0);
    test::assertEqual(m.get(1, 1).value(), 4.0);
}

void testDetermOne() {
    Matrix m(1);
    put(m, 0, 0, 5.0);
    test::assertAlmostEqual(m.determ().value(), 5.0);
}

void testDetermIdentity() {
    Matrix m(4);
    for (size_t i = 0; i < 4; ++i) put(m, i, i, 1.0);
    test::assertAlmostEqual(m.determ().value(), 1.0);
}

void testDetermZero() {
    Matrix m(3);
    test::assertAlmostEqual(m.determ().value(), 0.0);
}

void testDetermDiagonal() {
    Matrix m(4);
    put(m, 0, 0, 2.0); put(m, 1, 1, 3.0);
    put(m, 2, 2, 5.0); put(m, 3, 3, 7.0);
    test::assertAlmostEqual(m.determ().value(), 2.0 * 3.0 * 5.0 * 7.0);
}

void testDeterm2x2() {
    Matrix m(2);
    put(m, 0, 0, 1.0); put(m, 0, 1, 2.0);
    put(m, 1, 0, 3.0); put(m, 1, 1, 4.0);
    test::assertAlmostEqual(m.determ().value(), -2.0);     // 1·4 − 2·3
}

void testDeterm3x3() {
    // Expand along row 0: 6·(-38) − 1·(-14) + 36 = −178.
    Matrix m(3);
    put(m, 0, 0,  6.0); put(m, 0, 1,  1.0); put(m, 0, 2,  1.0);
    put(m, 1, 0,  4.0); put(m, 1, 1, -2.0); put(m, 1, 2,  5.0);
    put(m, 2, 0,  2.0); put(m, 2, 1,  8.0); put(m, 2, 2, -1.0);
    test::assertAlmostEqual(m.determ().value(), -178.0);
}

void testDeterm4x4UpperTri() {
    // Upper triangular: det = product of diagonal = 1·2·3·4 = 24.
    Matrix m(4);
    put(m, 0, 0, 1.0); put(m, 0, 1, 5.0); put(m, 0, 2, 9.0); put(m, 0, 3, 2.0);
                       put(m, 1, 1, 2.0); put(m, 1, 2, 6.0); put(m, 1, 3, 7.0);
                                          put(m, 2, 2, 3.0); put(m, 2, 3, 8.0);
                                                             put(m, 3, 3, 4.0);
    test::assertAlmostEqual(m.determ().value(), 24.0);
}

void testAdjugate2x2() {
    // adj([[a,b],[c,d]]) = [[d,-b],[-c,a]]
    Matrix m(2);
    put(m, 0, 0, 1.0); put(m, 0, 1, 2.0);
    put(m, 1, 0, 3.0); put(m, 1, 1, 4.0);
    Matrix adj = m.adjugate();
    test::assertAlmostEqual(adj.get(0, 0).value(),  4.0);
    test::assertAlmostEqual(adj.get(0, 1).value(), -2.0);
    test::assertAlmostEqual(adj.get(1, 0).value(), -3.0);
    test::assertAlmostEqual(adj.get(1, 1).value(),  1.0);
}

void testAdjugateIdentity() {
    Matrix m(3);
    for (size_t i = 0; i < 3; ++i) put(m, i, i, 1.0);
    Matrix adj = m.adjugate();
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            test::assertAlmostEqual(adj.get(i, j).value(), (i == j) ? 1.0 : 0.0);
}

// Verify M · adj(M) = |M| · I.
void testAdjugateRoundtrip3x3() {
    Matrix m(3);
    put(m, 0, 0,  6.0); put(m, 0, 1,  1.0); put(m, 0, 2,  1.0);
    put(m, 1, 0,  4.0); put(m, 1, 1, -2.0); put(m, 1, 2,  5.0);
    put(m, 2, 0,  2.0); put(m, 2, 1,  8.0); put(m, 2, 2, -1.0);
    const double det = m.determ().value();
    Matrix adj = m.adjugate();
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double acc = 0.0;
            for (size_t k = 0; k < 3; ++k)
                acc += m.get(i, k).value() * adj.get(k, j).value();
            test::assertAlmostEqual(acc, (i == j) ? det : 0.0, 1e-10);
        }
    }
}

void testDetermPropagatesVariance() {
    // 2x2 diagonal of independent uncertain entries — determinant gets non-zero unc.
    Matrix m(2);
    m.set(0, 0, VarDbl(3.0, 0.1));
    m.set(1, 1, VarDbl(5.0, 0.1));
    VarDbl d = m.determ();
    test::assertAlmostEqual(d.value(), 15.0);
    test::assertTrue(d.uncertainty() > 0.0, "uncertainty should propagate");
}

void testRandomMatrix() {
    const double dev = 1e-3;
    Matrix m = Matrix::randomMatrix(4, dev);
    const double expectedUnc = 256.0 * dev;
    for (size_t r = 0; r < 4; ++r) {
        for (size_t c = 0; c < 4; ++c) {
            test::assertAlmostEqual(m.get(r, c).uncertainty(), expectedUnc);
            const double v = m.get(r, c).value();
            test::assertTrue(v > -1024.0 && v < 1024.0,
                             "randomMatrix value out of expected range");
        }
    }
    VarDbl d = m.determ();
    (void) d;
}

// Single-pass matrix analysis equivalent to Python/runMatrixAnalysis.py:
// sweeps (size, noise) cells from size 4..9 across 20 noise levels
// (0 plus 1e-17..1e1), generates ~1000/size^2 random matrices plus one Hilbert
// per cell, and writes MatrixCondition_4_10.txt and AdjMatrix_4_10.txt under
// Cpp/Output. The conversion lives in Matrix::runMatrixAnalysis (Matrix.h);
// this entry point matches the Python driver's defaults.
void dumpMatrix() {
    const size_t minSize     = 4;
    const size_t maxSize     = 10;   // exclusive
    const size_t targetCells = 1000;
    const std::string outDir = "Output";

    std::cout << "Matrix analysis: sizes [" << minSize << ", " << maxSize
              << "), targetCells=" << targetCells
              << ", output=" << outDir << '\n';
    Matrix::runMatrixAnalysis(minSize, maxSize, targetCells,
                              Matrix::defaultNoises(), outDir);
    std::cout << "Matrix analysis complete\n";
}


int main(int argc, char* argv[]) {
    if (argc == 1) {
        testConstruction();
        testGetSet();
        testDetermOne();
        testDetermIdentity();
        testDetermZero();
        testDetermDiagonal();
        testDeterm2x2();
        testDeterm3x3();
        testDeterm4x4UpperTri();
        testAdjugate2x2();
        testAdjugateIdentity();
        testAdjugateRoundtrip3x3();
        testDetermPropagatesVariance();
        testRandomMatrix();
        std::cout << "All Matrix tests are successful";
    } else if ((argc == 2) && (std::string(argv[1]) == "Test"))
        dumpMatrix();
    return 0;
}
