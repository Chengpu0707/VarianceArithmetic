#include "IndexSin.h"
#include "Test.h"
#include "ulp.h"

#include <iostream>

using namespace var_dbl;

const IndexSin indexSin(3, false);

void test_neg_rem() // assume order==3, so size==8
{
    // when index is near 0
    test::assertEqual(0, 1 / 4);
    test::assertEqual(1, 1 % 4);

    test::assertEqual(0, -1 / 4);
    test::assertEqual(-1, -1 % 4);

    // when index is near +pi/2
    test::assertEqual(0, 3 / 4);    //0 % 3
    test::assertEqual(3, 3 % 4);
    test::assertEqual(1, 5 / 4);    //0 % 3
    test::assertEqual(1, 5 % 4);

    // when index is near +pi
    test::assertEqual(1, 7 / 4);    //0 % 1
    test::assertEqual(3, 7 % 4);
    test::assertEqual(2, 9 / 4);    //0 % -1
    test::assertEqual(1, 9 % 4);

    // when index is near -pi/2
    test::assertEqual(0, -3 / 4);   //0 % -3
    test::assertEqual(-3, -3 % 4);
    test::assertEqual(-1, -5 / 4);  //0 % -3
    test::assertEqual(-1, -5 % 4);

    // when index is near -pi
    test::assertEqual(-1, -7 / 4);  //0 % -1
    test::assertEqual(-3, -7 % 4);
    test::assertEqual(-2, -9 / 4);  //0 % 1
    test::assertEqual(-1, -9 % 4);
}

void test_get_sin_index()
{
    test::assertEqual(8, indexSin.size());

    test::assertEqual( 0, indexSin.get_sin_index(0));
    test::assertEqual( 1, indexSin.get_sin_index(1));
    test::assertEqual( 2, indexSin.get_sin_index(2));
    test::assertEqual( 3, indexSin.get_sin_index(3));
    test::assertEqual( 4, indexSin.get_sin_index(4));
    test::assertEqual( 3, indexSin.get_sin_index(5));
    test::assertEqual( 2, indexSin.get_sin_index(6));
    test::assertEqual( 1, indexSin.get_sin_index(7));
    test::assertEqual( 0, indexSin.get_sin_index(8));
    test::assertEqual(-1, indexSin.get_sin_index(9));
    test::assertEqual(-2, indexSin.get_sin_index(10));
    test::assertEqual(-3, indexSin.get_sin_index(11));
    test::assertEqual(-4, indexSin.get_sin_index(12));
    test::assertEqual(-3, indexSin.get_sin_index(13));
    test::assertEqual(-2, indexSin.get_sin_index(14));
    test::assertEqual(-1, indexSin.get_sin_index(15));
    test::assertEqual( 0, indexSin.get_sin_index(16));
    test::assertEqual( 1, indexSin.get_sin_index(17));
    test::assertEqual( 2, indexSin.get_sin_index(18));

    test::assertEqual(-1, indexSin.get_sin_index(-1));
    test::assertEqual(-2, indexSin.get_sin_index(-2));
    test::assertEqual(-3, indexSin.get_sin_index(-3));
    test::assertEqual(-4, indexSin.get_sin_index(-4));
    test::assertEqual(-3, indexSin.get_sin_index(-5));
    test::assertEqual(-2, indexSin.get_sin_index(-6));
    test::assertEqual(-1, indexSin.get_sin_index(-7));
    test::assertEqual( 0, indexSin.get_sin_index(-8));
    test::assertEqual( 1, indexSin.get_sin_index(-9));
    test::assertEqual( 2, indexSin.get_sin_index(-10));
    test::assertEqual( 3, indexSin.get_sin_index(-11));
    test::assertEqual( 4, indexSin.get_sin_index(-12));
    test::assertEqual( 3, indexSin.get_sin_index(-13));
    test::assertEqual( 2, indexSin.get_sin_index(-14));
    test::assertEqual( 1, indexSin.get_sin_index(-15));
    test::assertEqual( 0, indexSin.get_sin_index(-16));
    test::assertEqual(-1, indexSin.get_sin_index(-17));
    test::assertEqual(-2, indexSin.get_sin_index(-18));
}

void test_sin()
{
    test::assertEqual(8, indexSin.size());

    test::assertEqual(indexSin.sin(0).value(), 0);
    test::assertEqual(indexSin.sin(1).value(), std::sin(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(2).value(), std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.sin(3).value(), std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(4).value(), 1);
    test::assertEqual(indexSin.sin(5).value(), std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(6).value(), std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.sin(7).value(), std::sin(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(8).value(), 0);

    test::assertEqual(indexSin.sin(-1).value(), -std::sin(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(-2).value(), -std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.sin(-3).value(), -std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(-4).value(), -1);
    test::assertEqual(indexSin.sin(-5).value(), -std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(-6).value(), -std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.sin(-7).value(), -std::sin(IndexSin::PI*1/8));
    test::assertEqual(indexSin.sin(-8).value(), 0);
}

void test_sin_error()
{
    test::assertEquals(-std::sin(IndexSin::PI*1/8) - 9.992007221626409e-16, 
                        std::sin(IndexSin::PI*1001/8));
    test::assertEquals(ulp(std::sin(IndexSin::PI*1/8)), 5.551115123125783e-17);
}

void test_cos()
{
    test::assertEqual(indexSin.cos(0).value(), 1);
    test::assertEqual(indexSin.cos(1).value(), std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.cos(2).value(), std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.cos(3).value(), std::sin(IndexSin::PI*1/8));
    test::assertEqual(indexSin.cos(4).value(), 0);
    test::assertEqual(indexSin.cos(5).value(), -std::sin(IndexSin::PI*1/8));
    test::assertEqual(indexSin.cos(6).value(), -std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.cos(7).value(), -std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.cos(8).value(), -1);

    test::assertEqual(indexSin.cos(-1).value(), std::cos(IndexSin::PI*1/8));
    test::assertEqual(indexSin.cos(-2).value(), std::cos(IndexSin::PI*2/8));
    test::assertEqual(indexSin.cos(-3).value(), std::sin(IndexSin::PI*1/8));
}

int main()
{
    test_neg_rem();
    test_get_sin_index();
    test_sin();
    test_sin_error();
    test_cos();

    return 0;
}