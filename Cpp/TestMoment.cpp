/*
Unit tests for Moment.h NormalMoment: verifies bounding, leakage and computed
moment table for Gaussian distributions and writes a reference table to
Cpp/Output/ for cross-implementation comparison.
*/
#include "Moment.h"
#include "Test.h"
#include "ulp.h"

#include <fstream>


using namespace var_dbl;

int main()
{
    const NormalMoment mmt(5.0);
    test::assertEqual( mmt.bounding, 5);
    test::assertEqual( mmt.maxOrder(), 448);
    test::assertAlmostEqual( mmt.leakage, 5.7330e-07, 1e-10);
    // Normalized per Formula (2.2): ζ(0,κ)=1 by construction.
    test::assertAlmostEqual( mmt[0], 1.0, 1e-12);
    test::assertAlmostEqual( mmt[2] / 0.9999851327963293, 1., 1e-10);
    test::assertAlmostEqual( mmt[4] / 2.9995837182972154, 1., 1e-10);
    test::assertAlmostEqual( mmt[6] / 14.988626589191844, 1., 1e-10);
    test::assertAlmostEqual( mmt[8] / 104.68808606698785, 1., 1e-10);
    test::assertAlmostEqual( mmt[10] / 936.3852731689979, 1., 1e-10);

    std::ofstream ofs("../Cpp/Output/NormalMoment_5.0.txt");
    test::assertTrue(ofs.is_open(), "../Cpp/Output/NormalMoment_5.0.txt");
    ofs << "n\tMoment\n";
    for (int n = 0; n < mmt.maxOrder(); n += 2)
        ofs << n << "\t" << mmt[n] << "\n";

    std::cout << "All Moment tests are successful";
    return 0;
}
