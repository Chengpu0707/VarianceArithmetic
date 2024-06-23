#include "IndexSin.h"
#include "Stat.h"
#include "Test.h"
#include "ulp.h"

#include <fstream>
#include <iostream>
#include <numbers>

using namespace var_dbl;

const IndexSin indexSin(3);

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

void test_get_index()
{
    test::assertEqual(8, indexSin.size());

    test::assertEqual( 0, indexSin.get_index(0));
    test::assertEqual( 1, indexSin.get_index(1));
    test::assertEqual( 2, indexSin.get_index(2));
    test::assertEqual( 3, indexSin.get_index(3));
    test::assertEqual( 4, indexSin.get_index(4));
    test::assertEqual( 3, indexSin.get_index(5));
    test::assertEqual( 2, indexSin.get_index(6));
    test::assertEqual( 1, indexSin.get_index(7));
    test::assertEqual( 0, indexSin.get_index(8));
    test::assertEqual(-1, indexSin.get_index(9));
    test::assertEqual(-2, indexSin.get_index(10));
    test::assertEqual(-3, indexSin.get_index(11));
    test::assertEqual(-4, indexSin.get_index(12));
    test::assertEqual(-3, indexSin.get_index(13));
    test::assertEqual(-2, indexSin.get_index(14));
    test::assertEqual(-1, indexSin.get_index(15));
    test::assertEqual( 0, indexSin.get_index(16));
    test::assertEqual( 1, indexSin.get_index(17));
    test::assertEqual( 2, indexSin.get_index(18));

    test::assertEqual(-1, indexSin.get_index(-1));
    test::assertEqual(-2, indexSin.get_index(-2));
    test::assertEqual(-3, indexSin.get_index(-3));
    test::assertEqual(-4, indexSin.get_index(-4));
    test::assertEqual(-3, indexSin.get_index(-5));
    test::assertEqual(-2, indexSin.get_index(-6));
    test::assertEqual(-1, indexSin.get_index(-7));
    test::assertEqual( 0, indexSin.get_index(-8));
    test::assertEqual( 1, indexSin.get_index(-9));
    test::assertEqual( 2, indexSin.get_index(-10));
    test::assertEqual( 3, indexSin.get_index(-11));
    test::assertEqual( 4, indexSin.get_index(-12));
    test::assertEqual( 3, indexSin.get_index(-13));
    test::assertEqual( 2, indexSin.get_index(-14));
    test::assertEqual( 1, indexSin.get_index(-15));
    test::assertEqual( 0, indexSin.get_index(-16));
    test::assertEqual(-1, indexSin.get_index(-17));
    test::assertEqual(-2, indexSin.get_index(-18));
}

void test_sin()
{
    test::assertEqual(8, indexSin.size());

    test::assertAlmostEqual(indexSin.sin(0).value(), 0);
    test::assertAlmostEqual(indexSin.sin(1).value(), std::sin(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(2).value(), std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.sin(3).value(), std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(4).value(), 1);
    test::assertAlmostEqual(indexSin.sin(5).value(), std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(6).value(), std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.sin(7).value(), std::sin(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(8).value(), 0);

    test::assertAlmostEqual(indexSin.sin(-1).value(), -std::sin(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(-2).value(), -std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.sin(-3).value(), -std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(-4).value(), -1);
    test::assertAlmostEqual(indexSin.sin(-5).value(), -std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(-6).value(), -std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.sin(-7).value(), -std::sin(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.sin(-8).value(), 0);
}

void test_sin_error()
{
    test::assertAlmostEqual(-std::sin(std::numbers::pi*1/8) - 9.992007221626409e-16, 
                            std::sin(std::numbers::pi*1001/8));
    test::assertAlmostEqual(ulp(std::sin(std::numbers::pi*1/8)), 5.551115123125783e-17);
}

void test_cos()
{
    test::assertAlmostEqual(indexSin.cos(0).value(), 1);
    test::assertAlmostEqual(indexSin.cos(1).value(), std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.cos(2).value(), std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.cos(3).value(), std::sin(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.cos(4).value(), 0);
    test::assertAlmostEqual(indexSin.cos(5).value(), -std::sin(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.cos(6).value(), -std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.cos(7).value(), -std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.cos(8).value(), -1);

    test::assertAlmostEqual(indexSin.cos(-1).value(), std::cos(std::numbers::pi*1/8));
    test::assertAlmostEqual(indexSin.cos(-2).value(), std::cos(std::numbers::pi*2/8));
    test::assertAlmostEqual(indexSin.cos(-3).value(), std::sin(std::numbers::pi*1/8));
}

void dump_error()
{
    struct Hint {
        unsigned index;
        double rad;
        double first;
        double second;
        double ulp;
    };
    const unsigned max_order = 12;
    Stat<double, Hint> sin, cos, tan, lib;

    std::ostringstream oss;
    oss << "./Output/Lib_Error_" << max_order << ".txt";
    std::ofstream of(oss.str());
    of << "Order\tSize";
    for (auto subj: {"Sin", "Cos", "Tan", "Lib"}) {
        for (auto prop: {" Mean", " Std",  
                         " Max", " Max At Ulp", " Max At Index", " Max At Rad", " Max At First", " Max At Second", 
                         " Min", " Min At Ulp", " Min At Index", " Min At Rad", " Min At First", " Min At Second"})
            of << "\t" << subj << prop;
    }
    of << "\n";
    of.flush();
    for (unsigned order = 3; order < max_order; ++order) {
        const IndexSin indexSin(order);
        const unsigned size = 1 << (order * 2);
        for (unsigned i = 0; i < size; ++i) {
            const double isin = indexSin.sin(i).value();
            const double icos = indexSin.cos(i).value();
            const double itan = isin/icos;
            const double rad = std::numbers::pi * i / (1 << order);
            const double rsin = std::sin(rad);
            const double rcos = std::cos(rad);
            const double rtan = std::tan(rad);
            const double rlib = rsin / rcos;
            if ((rsin + isin) != 0) {
                const double dulp = (rsin - isin) / ulp(rsin + isin);
                sin.addAt(rsin - isin, {i, rad, isin, rsin, dulp});
            }
            if ((rcos + icos) != 0) {
                const double dulp = (rcos - icos) / ulp(rcos + icos);
                cos.addAt(rcos - icos, {i, rad, icos, rcos, dulp});
            }
            if (std::isfinite(itan) && std::isfinite(rtan)) {
                const double dulp = (itan - rtan) / ulp(itan + rtan);
                tan.addAt(dulp, {i, rad, itan, rtan, dulp});
            }
            if (std::isfinite(rlib) && std::isfinite(rtan)) {
                const double dulp = (rlib - rtan) / ulp(rlib + rtan);
                lib.addAt(dulp, {i, rad, rlib, rtan, dulp});
            }
        }
        of << order << "\t" << size;
        for (auto subj: {sin, cos, tan, lib}) {
            of << "\t" << subj.mean();
            of << "\t" << subj.std();
            of << "\t" << subj.max();
            Hint maxAt = subj.maxAt();
            of << "\t" << maxAt.ulp << "\t" << maxAt.index << "\t"  << maxAt.rad << "\t"  << maxAt.first << "\t"  << maxAt.second;
            of << "\t" << subj.min();
            Hint minAt = subj.minAt();
            of << "\t" << minAt.ulp << "\t" << minAt.index << "\t"  << minAt.rad << "\t"  << minAt.first << "\t"  << minAt.second;
        }
        of << "\n";
        of.flush();
    }
}


int main()
{
    test_neg_rem();
    test_get_index();
    test_sin();
    test_sin_error();
    test_cos();
    dump_error();
    std::cout << "All indexed sine tests are successful";
    return 0;
}