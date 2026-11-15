/*
Unit tests for Taylor::poly1d: verifies value/uncertainty propagation through
constant, linear, quadratic and higher-degree integer polynomials evaluated on
VarDbl inputs.
*/
#include "Taylor.h"
#include "Test.h"

using namespace var_dbl;


// Polynomial coefficients without C++11 initializer lists.
template <typename T> std::vector<T> coeff(T c0) {
    return std::vector<T>(1, c0);
}
template <typename T> std::vector<T> coeff(T c0, T c1) {
    std::vector<T> v = coeff(c0); v.push_back(c1); return v;
}
template <typename T> std::vector<T> coeff(T c0, T c1, T c2) {
    std::vector<T> v = coeff(c0, c1); v.push_back(c2); return v;
}
template <typename T> std::vector<T> coeff(T c0, T c1, T c2, T c3) {
    std::vector<T> v = coeff(c0, c1, c2); v.push_back(c3); return v;
}

void test_poly_0() {
    VarDbl res;
    try {
        res = Taylor::poly1d(VarDbl(), coeff(1));
        test::assertAlmostEqual(res.value(), 1., 0.);
        test::assertAlmostEqual(res.uncertainty(), 0., 0.);

        res = Taylor::poly1d(VarDbl(1), coeff(1));
        test::assertAlmostEqual(res.value(), 1., 0.);
        test::assertAlmostEqual(res.uncertainty(), 0., 0.);

        res = Taylor::poly1d(VarDbl(1), coeff(-2));
        test::assertAlmostEqual(res.value(), -2., 0.);
        test::assertAlmostEqual(res.uncertainty(), 0., 0.);

        res = Taylor::poly1d(VarDbl(2), coeff(1.));
        test::assertAlmostEqual(res.value(), 1., 0.);
        test::assertAlmostEqual(res.uncertainty(), 0., 0.);
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }
}


void test_poly_1() {
    VarDbl res;
    try {
        res = Taylor::poly1d(VarDbl(0, 1./8), coeff(0, 1));
        test::assertAlmostEqual(res.value(), 0., 0.);
        test::assertAlmostEqual(res.uncertainty(), 1./8, 1e-6);

        res = Taylor::poly1d(VarDbl(0, 1./8), coeff(1,1));
        test::assertAlmostEqual(res.value(), 1., 0.);
        test::assertAlmostEqual(res.uncertainty(), 1./8, 1e-6);

        res = Taylor::poly1d(VarDbl(-2, 1./8), coeff(1,1));
        test::assertAlmostEqual(res.value(), -1., 0.);
        test::assertAlmostEqual(res.uncertainty(), 1./8, 1e-6);
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }
}


void test_poly_2() {
    VarDbl res, res2;
    try {
        const double sValue[] = {0., 1., -1., 2., -2., 0.25, -0.25};
        for (size_t i = 0; i < sizeof(sValue) / sizeof(sValue[0]); ++i) {
            const double value = sValue[i];
            res = Taylor::poly1d(VarDbl(value, 0.5), coeff(0,0,1));
            res2 = VarDbl(value*value + 0.5*0.5, std::sqrt(4 *value*value *0.5*0.5 + 2 *0.5*0.5*0.5*0.5));
            test::assertAlmostEqual(res.value(), res2.value(), 1e-5);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 5e-5);

            res = Taylor::poly1d(VarDbl(value, 0.5), coeff(1,2,1));
            res2 = VarDbl((value + 1)*(value + 1) + 0.5*0.5, std::sqrt(4 *(value + 1)*(value + 1) *0.5*0.5 + 2 *0.5*0.5*0.5*0.5));
            test::assertAlmostEqual(res.value(), res2.value(), 1e-5);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 5e-5);

            res = Taylor::poly1d(VarDbl(value, 0.5), coeff(1,-2,1));
            res2 = VarDbl((value - 1)*(value - 1) + 0.5*0.5, std::sqrt(4 *(value - 1)*(value - 1) *0.5*0.5 + 2 *0.5*0.5*0.5*0.5));
            test::assertAlmostEqual(res.value(), res2.value(), 1e-5);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 5e-5);
        }
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }
}


void test_poly_3() {
    VarDbl res, res2;
    try {
        const double sValue[] = {1., -1., 2., -2., 0.25, -0.25};
        for (size_t i = 0; i < sizeof(sValue) / sizeof(sValue[0]); ++i) {
            const double value = sValue[i];
            res = Taylor::poly1d(VarDbl(value, 0.5), coeff(0,0,0,1));
            res2 = Taylor::pow(VarDbl(value, 0.5), 3);
            test::assertAlmostEqual(res.value(), res2.value(), 0.);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 0.);

            res = Taylor::poly1d(VarDbl(value, 0.5), coeff(1,3,3,1));
            res2 = Taylor::pow(VarDbl(1 + value, 0.5), 3);
            test::assertAlmostEqual(res.value(), res2.value(), 0.);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 0.);

            res = Taylor::poly1d(VarDbl(value, 0.5), coeff(1,-3,3,-1));
            res2 = Taylor::pow(VarDbl(1 - value, 0.5), 3);
            test::assertAlmostEqual(res.value(), res2.value(), 0.);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 0.);
        }
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }
}


int main() {
    test_poly_0();
    test_poly_1();
    test_poly_2();
    test_poly_3();

    std::cout << "All polynominal tests are successful";
    return 0;
}