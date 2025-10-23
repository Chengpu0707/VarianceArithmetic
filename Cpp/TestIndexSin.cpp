#include "IndexSin.h"
#include "Stat.h"
#include "Test.h"
#include "ulp.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numbers>

using namespace var_dbl;

const double q1 = std::sin(std::numbers::pi*1/8);
const double q2 = std::sin(std::numbers::pi*2/8);
const double q3 = std::cos(std::numbers::pi*1/8);


void test_int_size()
    // the size of index
{
    test::assertEqual(sizeof(int), 4);
    test::assertEqual(sizeof(long), 4);
    test::assertEqual(sizeof(long long), 8);
    
    test::assertEqual(sizeof(unsigned), 4);
    test::assertEqual(sizeof(unsigned long), 4);
    test::assertEqual(sizeof(unsigned long long), 8);

    test::assertEqual(sizeof(double), 8);
    test::assertEqual(sizeof(long double), 16);
}

void test_floating_size()
{
    test::assertEqual(sizeof(std::numbers::pi), 8);
    const long double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164;
    const long double diff = std::sin(PI/4) - (long double) std::sin(std::numbers::pi/4);
    test::assertEqual((double) diff, 4.10370522285763428272e-17);
    test::assertEqual(sizeof(std::numbers::pi), 8);
}

void test_neg_rem()
    // negative of positive reminder
{
    div_t res;
    res = std::div(1, 8);
    test::assertEqual(res.quot, 0);
    test::assertEqual(res.rem, 1);

    res = std::div(5, 8);
    test::assertEqual(res.quot, 0);
    test::assertEqual(res.rem, 5);

    res = std::div(9, 8);
    test::assertEqual(res.quot, 1);
    test::assertEqual(res.rem, 1);

    res = std::div(13, 8);
    test::assertEqual(res.quot, 1);
    test::assertEqual(res.rem, 5);

    res = std::div(-1, 8);
    test::assertEqual(res.quot, 0);
    test::assertEqual(res.rem, -1);

    res = std::div(-5, 8);
    test::assertEqual(res.quot, 0);
    test::assertEqual(res.rem, -5);

    res = std::div(-9, 8);
    test::assertEqual(res.quot, -1);
    test::assertEqual(res.rem, -1);

    res = std::div(-13, 8);
    test::assertEqual(res.quot, -1);
    test::assertEqual(res.rem, -5);

}


void test_Full_Fixed_index()
{
    IndexSin idxSin(IndexSin::SinSource::Full);
    
    test::assertEqual(idxSin.getIndex(1, 3), 1);
    test::assertEqual(idxSin.getIndex(4, 3), 4);
    test::assertEqual(idxSin.getIndex(5, 3), 5);
    test::assertEqual(idxSin.getIndex(9, 3), -1);
    test::assertEqual(idxSin.getIndex(12, 3), -4);
    test::assertEqual(idxSin.getIndex(13, 3), -5);
    test::assertEqual(idxSin.getIndex(17, 3), 1);
    
    test::assertEqual(idxSin.getIndex(-1, 3), -1);
    test::assertEqual(idxSin.getIndex(-4, 3), -4);
    test::assertEqual(idxSin.getIndex(-5, 3), -5);
    test::assertEqual(idxSin.getIndex(-9, 3), 1);
    test::assertEqual(idxSin.getIndex(-12, 3), 4);
    test::assertEqual(idxSin.getIndex(-13, 3), 5);
    test::assertEqual(idxSin.getIndex(-17, 3), -1);
}


void test_Quart_index()
{
    IndexSin idxSin(IndexSin::SinSource::Quart);

    test::assertEqual(idxSin.getIndex(1, 3),    1);
    test::assertEqual(idxSin.getIndex(5, 3),    3);
    test::assertEqual(idxSin.getIndex(9, 3),   -1);
    test::assertEqual(idxSin.getIndex(13, 3),  -3);
    test::assertEqual(idxSin.getIndex(17, 3),   1);
       
    test::assertEqual(idxSin.getIndex(-1, 3),  -1);
    test::assertEqual(idxSin.getIndex(-5, 3),  -3);
    test::assertEqual(idxSin.getIndex(-9, 3),   1);
    test::assertEqual(idxSin.getIndex(-13, 3),  3);
    test::assertEqual(idxSin.getIndex(-17, 3), -1);
}



void test_Quart_sin()
{
    const IndexSin indexSin;

    test::assertAlmostEqual(indexSin.sin(0, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(1, 3).value(), q1);
    test::assertAlmostEqual(indexSin.sin(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(3, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(4, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.sin(5, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(6, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(7, 3).value(), q1);
    test::assertAlmostEqual(indexSin.sin(8, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(9, 3).value(), -q1);

    test::assertAlmostEqual(indexSin.sin(-1, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.sin(-2, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-3, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-4, 3).value(), -1.);
    test::assertAlmostEqual(indexSin.sin(-5, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-7, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.sin(-8, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(-9, 3).value(), q1);
}

void test_Quart_cos()
{
    const IndexSin indexSin;

    test::assertAlmostEqual(indexSin.cos(0, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.cos(1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(3, 3).value(), q1);
    test::assertAlmostEqual(indexSin.cos(4, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.cos(5, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.cos(6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.cos(7, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.cos(8, 3).value(), -1.);

    test::assertAlmostEqual(indexSin.cos(-1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(-2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(-3, 3).value(), q1);
}

void test_Full_sin()
{
    const IndexSin indexSin(IndexSin::SinSource::Full);

    test::assertAlmostEqual(indexSin.sin(0, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(1, 3).value(), q1);
    test::assertAlmostEqual(indexSin.sin(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(3, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(4, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.sin(5, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(6, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(7, 3).value(), q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.sin(8, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(9, 3).value(), -q1);

    test::assertAlmostEqual(indexSin.sin(-1, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.sin(-2, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-3, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-4, 3).value(), -1.);
    test::assertAlmostEqual(indexSin.sin(-5, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-7, 3).value(), -q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.sin(-8, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(-9, 3).value(), q1);
}

void test_Full_cos()
{
    const IndexSin indexSin(IndexSin::SinSource::Full);

    IndexSin idxSin(IndexSin::SinSource::Full);
    test::assertAlmostEqual(indexSin.cos(0, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.cos(1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(3, 3).value(), q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.cos(4, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.cos(5, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.cos(6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.cos(7, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.cos(8, 3).value(), -1.);

    test::assertAlmostEqual(indexSin.cos(-1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(-2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(-3, 3).value(), q1);
}

void dump_Quart_indexSin()
{
    const IndexSin indexSin(IndexSin::SinSource::Quart);
    test::assertTrue(indexSin.dump(10, "./Output/IndexSin_Quart_10.txt"));
    const IndexSin readback(IndexSin::SinSource::Quart, "./Output");
    const unsigned order = 6;
    for (size_t n = 0; n <= (1 << order); ++n) {
        test::assertAlmostEqual(indexSin.sin(n, order).value(), readback.sin(n, order).value());
        test::assertAlmostEqual(indexSin.sin(n, order).uncertainty(), readback.sin(n, order).uncertainty());
        test::assertAlmostEqual(indexSin.cos(n, order).value(), readback.cos(n, order).value());
        test::assertAlmostEqual(indexSin.cos(n, order).uncertainty(), readback.cos(n, order).uncertainty());
    }
    const IndexSin python(IndexSin::SinSource::Quart, "../Python/Output");
    for (size_t n = 0; n <= (1 << order); ++n) {
        test::assertAlmostEqual(indexSin.sin(n, order).value(), python.sin(n, order).value());
        test::assertAlmostEqual(indexSin.sin(n, order).uncertainty(), python.sin(n, order).uncertainty());
        test::assertAlmostEqual(indexSin.cos(n, order).value(), python.cos(n, order).value());
        test::assertAlmostEqual(indexSin.cos(n, order).uncertainty(), python.cos(n, order).uncertainty());
    }
}

void dump_Full_indexSin()
{
    const IndexSin indexSin(IndexSin::SinSource::Full);
    test::assertTrue(indexSin.dump(10, "./Output/IndexSin_Full_10.txt"));
    const IndexSin readback(IndexSin::SinSource::Full, "./Output");
    const unsigned order = 6;
    for (size_t n = 0; n <= (1 << order); ++n) {
        test::assertAlmostEqual(indexSin.sin(n, order).value(), readback.sin(n, order).value());
        test::assertAlmostEqual(indexSin.cos(n, order).uncertainty(), readback.cos(n, order).uncertainty());
    }
    const IndexSin python(IndexSin::SinSource::Full, "../Python/Output");
    for (size_t n = 0; n <= (1 << order); ++n) {
        test::assertAlmostEqual(indexSin.sin(n, order).value(), python.sin(n, order).value());
        test::assertAlmostEqual(indexSin.cos(n, order).uncertainty(), python.cos(n, order).uncertainty());
    }
}

void test_Lib_sin()
{
    const IndexSin indexSin(IndexSin::SinSource::Lib);

    test::assertAlmostEqual(indexSin.sin(0, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(1, 3).value(), q1);
    test::assertAlmostEqual(indexSin.sin(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(3, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(4, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.sin(5, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(6, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(7, 3).value(), q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.sin(8, 3).value(), 0., ulp(1.));
    test::assertAlmostEqual(indexSin.sin(9, 3).value(), -q1, ulp(q1)*2);

    test::assertAlmostEqual(indexSin.sin(-1, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.sin(-2, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-3, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-4, 3).value(), -1.);
    test::assertAlmostEqual(indexSin.sin(-5, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-7, 3).value(), -q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.sin(-8, 3).value(), 0., ulp(1.));
    test::assertAlmostEqual(indexSin.sin(-9, 3).value(), q1, ulp(q1)*2);

    test::assertAlmostEqual(indexSin.sin(9, 3).value(), -q1, ulp(q1)*2);    
    test::assertAlmostEqual(indexSin.sin(97, 3).value(), q1, ulp(q1)*63);    
    test::assertAlmostEqual(indexSin.sin(1001, 3).value(), -q1, ulp(q1)*18);    
    test::assertAlmostEqual(indexSin.sin(10001, 3).value(), q1, ulp(q1)*4108);    
}

void test_Lib_cos()
{
    const IndexSin indexSin(IndexSin::SinSource::Lib);

    test::assertAlmostEqual(indexSin.cos(0, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.cos(1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(3, 3).value(), q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.cos(4, 3).value(), 0., ulp(1.));
    test::assertAlmostEqual(indexSin.cos(5, 3).value(), -q1, ulp(q1)*2);
    test::assertAlmostEqual(indexSin.cos(6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.cos(7, 3).value(), -q3, ulp(q3)*2);
    test::assertAlmostEqual(indexSin.cos(8, 3).value(), -1.);

    test::assertAlmostEqual(indexSin.cos(-1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(-2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(-3, 3).value(), q1);
}


void test_Prec_index()
{
    IndexSin idxSin(IndexSin::SinSource::Prec);

    test::assertEqual(idxSin.getIndex(1, 3),    1);
    test::assertEqual(idxSin.getIndex(5, 3),    3);
    test::assertEqual(idxSin.getIndex(9, 3),   -1);
    test::assertEqual(idxSin.getIndex(13, 3),  -3);
    test::assertEqual(idxSin.getIndex(17, 3),   1);
       
    test::assertEqual(idxSin.getIndex(-1, 3),  -1);
    test::assertEqual(idxSin.getIndex(-5, 3),  -3);
    test::assertEqual(idxSin.getIndex(-9, 3),   1);
    test::assertEqual(idxSin.getIndex(-13, 3),  3);
    test::assertEqual(idxSin.getIndex(-17, 3), -1);
}

void test_Prec_sin()
{
    const IndexSin indexSin(IndexSin::SinSource::Prec);

    test::assertAlmostEqual(indexSin.sin(0, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(1, 3).value(), q1);
    test::assertAlmostEqual(indexSin.sin(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(3, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(4, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.sin(5, 3).value(), q3);
    test::assertAlmostEqual(indexSin.sin(6, 3).value(), q2);
    test::assertAlmostEqual(indexSin.sin(7, 3).value(), q1);
    test::assertAlmostEqual(indexSin.sin(8, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(9, 3).value(), -q1);

    test::assertAlmostEqual(indexSin.sin(-1, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.sin(-2, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-3, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-4, 3).value(), -1.);
    test::assertAlmostEqual(indexSin.sin(-5, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.sin(-6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.sin(-7, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.sin(-8, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.sin(-9, 3).value(), q1);
}

void test_Prec_cos()
{
    const IndexSin indexSin(IndexSin::SinSource::Prec);

    test::assertAlmostEqual(indexSin.cos(0, 3).value(), 1.);
    test::assertAlmostEqual(indexSin.cos(1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(3, 3).value(), q1);
    test::assertAlmostEqual(indexSin.cos(4, 3).value(), 0.);
    test::assertAlmostEqual(indexSin.cos(5, 3).value(), -q1);
    test::assertAlmostEqual(indexSin.cos(6, 3).value(), -q2);
    test::assertAlmostEqual(indexSin.cos(7, 3).value(), -q3);
    test::assertAlmostEqual(indexSin.cos(8, 3).value(), -1.);

    test::assertAlmostEqual(indexSin.cos(-1, 3).value(), q3);
    test::assertAlmostEqual(indexSin.cos(-2, 3).value(), q2);
    test::assertAlmostEqual(indexSin.cos(-3, 3).value(), q1);
}

void dump_prec_indexSin()
{
    const IndexSin indexSin(IndexSin::SinSource::Prec);
    indexSin.dump(IndexSin::MAX_ORDER, "./Output/IndexSin_Prec.txt");
}


int main()
{
    test_int_size();
    test_floating_size();
    test_neg_rem();
    
    test_Prec_index();
    test_Prec_sin();
    test_Prec_cos();
    dump_prec_indexSin();

    test_Quart_index();
    test_Quart_sin();
    test_Quart_cos();
    dump_Quart_indexSin();

    test_Full_Fixed_index();
    test_Full_sin();
    test_Full_cos();
    dump_Full_indexSin();

    test_Lib_sin();
    test_Lib_cos();

    std::cout << "All indexed sine tests are successful";
    return 0;
}