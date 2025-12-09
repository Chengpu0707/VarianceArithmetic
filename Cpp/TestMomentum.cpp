#include "Momentum.h"
#include "Test.h"
#include "ulp.h"

#include <fstream>


using namespace var_dbl;

int main() 
{
    const NormalMomentum mmt(5.0);
    test::assertEqual( mmt.bounding, 5);
    test::assertEqual( mmt.maxOrder(), 448);
    test::assertAlmostEqual( mmt.leakage, 5.7330e-07, 1e-10);
    test::assertAlmostEqual( mmt[0] / 0.999999426696856, 1., 1e-10);
    test::assertAlmostEqual( mmt[2] / 0.999984559501709, 1., 1e-10);
    test::assertAlmostEqual( mmt[4] / 2.99958199862644, 1., 1e-10);
    test::assertAlmostEqual( mmt[6] / 14.9886179961651, 1., 1e-10);
    test::assertAlmostEqual( mmt[8] / 104.688026048979, 1., 1e-10);
    test::assertAlmostEqual( mmt[10] / 936.384736336377, 1., 1e-10);

    std::ofstream ofs("../Cpp/Output/NormalMomentum_5.0.txt");
    test::assertTrue(ofs.is_open(), "../Cpp/Output/NormalMomentum_5.0.txt");
    ofs << "n\tMomentum\n";
    for (int n = 0; n < mmt.maxOrder(); n += 2)
        ofs << n << "\t" << mmt[n] << "\n";

    std::cout << "All Momentum tests are successful";
    return 0;
}