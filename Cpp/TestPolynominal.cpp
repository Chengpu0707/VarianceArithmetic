#include "Taylor.h"
#include "Test.h"

using namespace var_dbl;

void test_poly_0() {
    VarDbl res;
    try {
        res = Taylor::poly1d(VarDbl(), std::vector<int>{1});
        test::assertAlmostEqual(res.value(), 1, 0);
        test::assertAlmostEqual(res.uncertainty(), 0, 0);

        res = Taylor::poly1d(VarDbl(1), std::vector<int>{1});
        test::assertAlmostEqual(res.value(), 1, 1);
        test::assertAlmostEqual(res.uncertainty(), 0, 0);

        res = Taylor::poly1d(VarDbl(1), std::vector<int>{-2});
        test::assertAlmostEqual(res.value(), -2, 0);
        test::assertAlmostEqual(res.uncertainty(), 0, 0);

        res = Taylor::poly1d(VarDbl(2), std::vector<double>{1.0});
        test::assertAlmostEqual(res.value(), 1, 0);
        test::assertAlmostEqual(res.uncertainty(), 0, 0);
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }
}


void test_poly_1() {
    VarDbl res;
    try {
        res = Taylor::poly1d(VarDbl(0, 1/8), std::vector<int>{0, 1});
        test::assertAlmostEqual(res.value(), 0, 0);
        test::assertAlmostEqual(res.uncertainty(), 1/8, 1e-6);

        res = Taylor::poly1d(VarDbl(0, 1/8), std::vector<int>{1,1});
        test::assertAlmostEqual(res.value(), 1, 0);
        test::assertAlmostEqual(res.uncertainty(), 1/8, 1e-6);

        res = Taylor::poly1d(VarDbl(-2, 1/8), std::vector<int>{1,1});
        test::assertAlmostEqual(res.value(), -1, 0);
        test::assertAlmostEqual(res.uncertainty(), 1/8, 1e-6);
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }
}


void test_poly_2() {
    VarDbl res, res2;
    try {
        for (double value: {0., 1., -1., 2., -2., 0.25, -0.25}) {
            res = Taylor::poly1d(VarDbl(value, 0.5), std::vector<int>{0,0,1});
            res2 = VarDbl(value*value + 0.5*0.5, std::sqrt(4 *value*value *0.5*0.5 + 2 *0.5*0.5*0.5*0.5));
            test::assertAlmostEqual(res.value(), res2.value(), 1e-5);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 5e-5);

            res = Taylor::poly1d(VarDbl(value, 0.5), std::vector<int>{1,2,1});
            res2 = VarDbl((value + 1)*(value + 1) + 0.5*0.5, std::sqrt(4 *(value + 1)*(value + 1) *0.5*0.5 + 2 *0.5*0.5*0.5*0.5));
            test::assertAlmostEqual(res.value(), res2.value(), 1e-5);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 5e-5);

            res = Taylor::poly1d(VarDbl(value, 0.5), std::vector<int>{1,-2,1});
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
        for (double value: {1., -1., 2., -2., 0.25, -0.25}) {
            res = Taylor::poly1d(VarDbl(value, 0.5), std::vector<int>{0,0,0,1});
            res2 = Taylor::pow(VarDbl(value, 0.5), 3);
            test::assertAlmostEqual(res.value(), res2.value(), 0);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 0);

            res = Taylor::poly1d(VarDbl(value, 0.5), std::vector<int>{1,3,3,1});
            res2 = Taylor::pow(VarDbl(1 + value, 0.5), 3);
            test::assertAlmostEqual(res.value(), res2.value(), 0);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 0);

            res = Taylor::poly1d(VarDbl(value, 0.5), std::vector<int>{1,-3,3,-1});
            res2 = Taylor::pow(VarDbl(1 - value, 0.5), 3);
            test::assertAlmostEqual(res.value(), res2.value(), 0);
            test::assertAlmostEqual(res.uncertainty(), res2.uncertainty(), 0);
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