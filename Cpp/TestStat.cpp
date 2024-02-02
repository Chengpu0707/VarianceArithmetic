#include "Stat.h"
#include "Test.h"

#include <iostream>
#include <utility>

using namespace var_dbl;

void testStat()
{
    std::vector sData{1,3,2};
    Stat stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 3);
    test::assertEquals(stat.mean, 2);
    test::assertEquals(stat.stddev, 1);

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 2);
    test::assertEquals(stat.mean, 2);
    test::assertEquals(stat.stddev, std::sqrt(2));

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 1);
    test::assertEquals(stat.mean, 1);
    test::assertEquals(stat.stddev, std::nan("2"));

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    test::assertEquals(stat.count, 0);
    test::assertEquals(stat.mean, std::nan("1"));
    test::assertEquals(stat.stddev, std::nan("2"));
}

void testHisto()
{
    std::vector<int> sData;
    for (int i = -10; i <= 10; ++i) 
        sData.push_back(i);
    Histo histo = calcHisto(sData.begin(), sData.end(), 1.5, 2);
    test::assertEqual(histo.stat.count, 21);
    test::assertEquals(histo.stat.mean, 0);
    test::assertEquals(histo.stat.stddev, 10/std::sqrt(3)*1.0677078252031309);
    test::assertEquals(histo.sHisto.size(), 2*3);
    test::assertEqual(histo.less, 1);
    test::assertEqual(histo.more, 1);
    double sCenter[] = {-1.25f, -0.75f, -0.25f, 0.25f, 0.75f, 1.25f};
    for (int i = -3; i < 3; ++i) {
        test::assertEquals(histo.sHisto[i + 3].first, sCenter[i + 3]);
        test::assertEquals(histo.sHisto[i + 3].second, ((i == 0)? 4.0f : 3.0f)/19);
    }
}

int main() 
{
    testStat();
    testHisto();
    std::cout << "All Stat tests are successful";
}