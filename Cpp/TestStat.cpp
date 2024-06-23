#include "Stat.h"
#include "Test.h"

using namespace var_dbl;

void testStat()
{
    const std::vector<int> sSample{1, 2, 3};
    Stat<int> stat1(sSample.begin(), sSample.end());
    test::assertEqual(3, stat1.count());
    test::assertEqual(1, stat1.min());
    test::assertEqual(3, stat1.max());
    test::assertAlmostEqual(2, stat1.mean());
    test::assertAlmostEqual(1, stat1.std());

    Stat<int> stat2{1, 2, 3};
    test::assertEqual(3, stat2.count());
    test::assertEqual(1, stat2.min());
    test::assertEqual(3, stat2.max());
    test::assertAlmostEqual(2, stat2.mean());
    test::assertAlmostEqual(1, stat2.std());

    Stat<double> stat3;
    test::assertEqual(0, stat3.count());
    test::assertTrue(stat3.min() > stat3.max());
    test::assertFalse(std::isfinite(stat3.mean()));
    test::assertFalse(std::isfinite(stat3.std()));

    test::assertEqual(1, stat3.addAt(1, 0));
    test::assertEqual(1, stat3.count());
    test::assertEqual(1, stat3.min());
    test::assertEqual(1, stat3.max());
    test::assertEqual(0, stat3.minAt());
    test::assertEqual(0, stat3.maxAt());
    test::assertAlmostEqual(1, stat3.mean());
    test::assertAlmostEqual(0, stat3.std());

    test::assertEqual(2, stat3.add(3));
    test::assertEqual(2, stat3.count());
    test::assertEqual(1, stat3.min());
    test::assertEqual(3, stat3.max());
    test::assertEqual(0, stat3.minAt());
    test::assertEqual(0, stat3.maxAt());
    test::assertAlmostEqual(2, stat3.mean());
    test::assertAlmostEqual(std::sqrt(2), stat3.std(), 1e-16);

    test::assertEqual(3, stat3.addAt(2, 1));
    test::assertEqual(3, stat3.count());
    test::assertEqual(1, stat3.min());
    test::assertEqual(3, stat3.max());
    test::assertEqual(0, stat3.minAt());
    test::assertEqual(0, stat3.maxAt());
    test::assertAlmostEqual(2, stat3.mean());
    test::assertAlmostEqual(1, stat3.std());

    test::assertEqual(6, stat3.add({1,2,3}));
    test::assertEqual(6, stat3.count());
    test::assertEqual(1, stat3.min());
    test::assertEqual(3, stat3.max());
    test::assertEqual(0, stat3.minAt());
    test::assertEqual(0, stat3.maxAt());
    test::assertAlmostEqual(2, stat3.mean());
    test::assertAlmostEqual(std::sqrt(0.8), stat3.std(), 3e-17);
}


void testHistogram()
{
    Histogram histo;
    test::assertEqual(histo.range, 3);
    test::assertEqual(histo.divids, 5);
    std::vector<double> sData;
    for (int i = -16; i <= 16; ++i) {
        test::assertEqual(histo.add(i/5.), i + 17);
        sData.push_back(i/5.);
    }
    test::assertEquals(histo.histogram(), std::vector<unsigned>(31, 1));
    test::assertEqual(histo.lowers(), 1);
    test::assertEqual(histo.uppers(), 1);

    histo.add(sData.begin(), sData.end());
    test::assertEquals(histo.histogram(), std::vector<unsigned>(31, 2));
    test::assertEqual(histo.lowers(), 2);
    test::assertEqual(histo.uppers(), 2);

    test::assertAlmostEqual(histo.mean(), 0);
    test::assertAlmostEqual(histo.std(), 16./5 /std::sqrt(3), 0.075);
        // 1.919 vs 1.849
}


void testWhite() 
{
    Random rand(0, 1);
    Histogram histo(std::sqrt(3), 2);
    for (unsigned i = 1; i <= 1000; ++i)
        test::assertEqual(histo.add(rand.white()), i);
    test::assertEqual(histo.lowers(), 0);
    test::assertEqual(histo.uppers(), 0);
    test::assertAlmostEqual(histo.mean(), 0, 0.07);
    test::assertAlmostEqual(histo.std(), 1, 0.03);
    const std::vector<unsigned> histogram = histo.histogram();
    test::assertEqual(histogram.size(), 7);
    for (unsigned i = 1; i < 4; ++i)
        test::assertAlmostEqual(histogram[i], 1000/7., 20);
    test::assertEqual(histo.header(), "\t-1.5\t-1\t-0.5\t0\t0.5\t1\t1.5");
    std::istringstream iss(histo.formatted());
    double pct;
    for (unsigned i = 0; i < 7; ++i) {
        iss >> pct;
        test::assertAlmostEqual(pct, 1./7, 0.05);
    }
}


void testGaussian() 
{
    Random rand(0, 1);
    Histogram sHist[]{Histogram(1, 10), Histogram(2, 10), Histogram(3, 10)};
    for (unsigned i = 1; i <= 10000; ++i) {
        const double val = rand.gauss();
        for (int j = 0; j < 3; ++j)
            test::assertEqual(sHist[j].add(val), i);
    }
    for (int j = 0; j < 3; ++j) {
        test::assertAlmostEqual(sHist[j].mean(), 0, 0.03);
        test::assertAlmostEqual(sHist[j].std(), 1, 0.02);
    }
    test::assertAlmostEqual(sHist[0].lowers() + sHist[0].uppers(), 3200, 1000);
    test::assertAlmostEqual(sHist[1].lowers() + sHist[1].uppers(), 500, 200);
    test::assertAlmostEqual(sHist[2].lowers() + sHist[2].uppers(), 30, 20);
}


int main()
{
    testStat();
    testHistogram();
    testWhite();
    testGaussian();
    std::cout << "All Stat tests are successful";
    return 0;
}