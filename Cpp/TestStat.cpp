#include "Stat.h"
#include "Test.h"

#include <iostream>
#include <utility>

using namespace var_dbl;

void testStat()
{
    std::vector sData{1,3,2};
    Stat stat = calcStat(sData.begin(), sData.end());
    Test::assertEquals(stat.count, 3);
    Test::assertEquals(stat.mean, 2);
    Test::assertEquals(stat.stddev, 1);

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    Test::assertEquals(stat.count, 2);
    Test::assertEquals(stat.mean, 2);
    Test::assertEquals(stat.stddev, std::sqrt(2));

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    Test::assertEquals(stat.count, 1);
    Test::assertEquals(stat.mean, 1);
    Test::assertEquals(stat.stddev, std::nan("2"));

    sData.pop_back();
    stat = calcStat(sData.begin(), sData.end());
    Test::assertEquals(stat.count, 0);
    Test::assertEquals(stat.mean, std::nan("1"));
    Test::assertEquals(stat.stddev, std::nan("2"));
}

void testHisto()
{
    std::vector<int> sData;
    for (int i = -10; i <= 10; ++i) 
        sData.push_back(i);
    Histo histo = calcHisto(sData.begin(), sData.end(), 1.5, 2);
    Test::assertEqual(histo.stat.count, 21);
    Test::assertEquals(histo.stat.mean, 0);
    Test::assertEquals(histo.stat.stddev, 10/std::sqrt(3)*1.0677078252031309);
    Test::assertEquals(histo.sHisto.size(), 2*3);
    Test::assertEqual(histo.less, 1);
    Test::assertEqual(histo.more, 1);
    double sCenter[] = {-1.25f, -0.75f, -0.25f, 0.25f, 0.75f, 1.25f};
    for (int i = -3; i < 3; ++i) {
        Test::assertEquals(histo.sHisto[i + 3].first, sCenter[i + 3]);
        Test::assertEquals(histo.sHisto[i + 3].second, ((i == 0)? 4.0f : 3.0f)/19);
    }
}

int main() 
{
    testStat();
    testHisto();
    std::cout << "All Stat tests are successful";
}