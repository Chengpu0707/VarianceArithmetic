#include "Momentum.h"
#include "Test.h"
#include "ulp.h"

#include <fstream>


using namespace var_dbl;

int main() 
{
    const auto _momentum = NormalMomentum();
    test::assertEqual( _momentum.size(), 442);
    test::assertEqual(  1 * 0.99999942673466358, _momentum[0], "0");
    test::assertEqual(  1 * 0.99998456037129324, _momentum[2], "2");
    test::assertEqual(  3 * 0.99986067282526425, _momentum[4], "4");
    test::assertEqual( 15 * 0.99924122967608053, _momentum[6], "6");
    test::assertEqual(105 * 0.99702891516115077, _momentum[8], "8");
    test::assertEqual(945 * 0.99088355330473754, _momentum[10], "10"); 

    test::assertEqual(0, _momentum[1], "1");
    test::assertEqual(0, _momentum[3], "3");
    test::assertEqual(0, _momentum[5], "5");
    test::assertEqual(0, _momentum[7], "7");
    test::assertEqual(0, _momentum[9], "9");

    std::ifstream ifs("../Python/Output/NormalMomentum_5.txt");
    test::assertTrue(ifs.is_open(), "../Python/Output/NormalMomentum_5.txt");
    std::string line;
    std::getline(ifs, line);
    test::assertEqual(line, "n\tMomentum\t!!Diff\tSigma=5.0", line);

    std::ofstream ofs("./Output/NormalMomentum_5.txt");
    test::assertTrue(ofs.is_open(), "./Output/NormalMomentum_5.txt");
    ofs << "2n\tPython\tCpp\tError" << std::endl;

    size_t j;
    double val;
    for (size_t i =0; i < _momentum.size(); i += 2) {
        ifs >> j >> val;
        test::assertEqual(i, j);
        test::assertAlmostEqual(_momentum[i]/val, 1, 2e-2);
        ofs << i << '\t' << val << '\t' << _momentum[i] << '\t' << _momentum[i]/val - 1 << std::endl;
    }

    std::cout << "All Momentum tests are successful";
    return 0;
}