#include "Momentum.h"
#include "Test.h"
#include "ulp.h"

#include <fstream>


using namespace var_dbl;

int main() 
{
    const Momentum& mmt = IDEAL_MOMENTUM;
    test::assertAlmostEqual( mmt.bounding, 5);
    test::assertAlmostEqual( mmt.leakage, 5.7330e-07, 1e-10);
    test::assertEqual( mmt.maxOrder(), 448);

    std::ifstream ifs("../Python/NormalMomentum_5.0.txt");
    test::assertTrue(ifs.is_open(), "../Python/NormalMomentum_5.0.txt");
    std::string line;
    std::getline(ifs, line);
    test::assertEqual(line, "n\tMomentum\tBounding:\t5.0", line);
    size_t j;
    double val;
    for (size_t i = 0; i <= 24; i += 2) {
        ifs >> j >> val;
        test::assertEqual(i, j);
        test::assertAlmostEqual(mmt[i]/val, 1, 1e-10);
        test::assertAlmostEqual(mmt[i + 1], 0, 1e-10);
    }

    std::ofstream ofs("../Cpp/Output/NormalMomentum_5.0.txt");
    test::assertTrue(ofs.is_open(), "../Python/Output/NormalMomentum_5.0.txt");
    ofs << "n\tMomentum\n";
    for (int n = 0; n < mmt.maxOrder(); n += 2)
        ofs << n << "\t" << mmt[n] << "\n";

    std::cout << "All Momentum tests are successful";
    return 0;
}