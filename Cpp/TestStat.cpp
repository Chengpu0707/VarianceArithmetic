#include "Stat.h"
#include "Test.h"

#include <iostream>

using namespace var_dbl;

void testStat()
{
    std::vector sData{1,3,2};
    Stat stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 3);
    test::assertEquals(stat.min, 1);
    test::assertEquals(stat.max, 3);
    test::assertEquals(stat.mean, 2);
    test::assertEquals(stat.stddev, 1);

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 2);
    test::assertEquals(stat.min, 1);
    test::assertEquals(stat.max, 3);
    test::assertEquals(stat.mean, 2);
    test::assertEquals(stat.stddev, std::sqrt(2));

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 1);
    test::assertEquals(stat.min, 1);
    test::assertEquals(stat.max, 1);
    test::assertEquals(stat.mean, 1);
    test::assertEquals(stat.stddev, std::numeric_limits<double>::quiet_NaN());

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.min, std::numeric_limits<double>::max());
    test::assertEquals(stat.max, -std::numeric_limits<double>::max());
    test::assertEquals(stat.count, 0);
    test::assertEquals(stat.mean, std::numeric_limits<double>::quiet_NaN());
    test::assertEquals(stat.stddev, std::numeric_limits<double>::quiet_NaN());
}

void testHisto()
{
    float RANGE = 1.5;
    unsigned DIVIDES = 2;
    std::vector<int> sData;
    for (int i = -10; i <= 10; ++i) 
        sData.push_back(i);
    Histo histo;
    test::assertTrue(calcHisto(histo, sData.begin(), sData.end(), RANGE, DIVIDES));
    test::assertEqual(histo.stat.count, 21);
    test::assertEqual(histo.stat.min, -10);
    test::assertEqual(histo.stat.max, 10);
    test::assertEquals(histo.stat.mean, 0);
    test::assertEquals(histo.stat.stddev, 10/std::sqrt(3)*1.0747092630102337);
    test::assertEquals(histo.sHisto.size(), 6);
    test::assertEqual(histo.less, 1);
    test::assertEqual(histo.more, 1);
    double sCenter[] = {-1.25f, -0.75f, -0.25f, 0.25f, 0.75f, 1.25f};
    std::vector<float> sHistoCenter = calcHistoCenters(RANGE, DIVIDES);
    test::assertEqual(sHistoCenter.size(), 6);
    for (int i = -3; i < 3; ++i) {
        test::assertEquals(sHistoCenter[i + 3], sCenter[i + 3]);
        test::assertEquals(histo.sHisto[i + 3], ((i == 0)? 4.0f : 3.0f)/19);
    }
}

int main() 
{
    testStat();
    testHisto();
    std::cout << "All Stat tests are successful";
}