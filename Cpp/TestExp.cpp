#include "TestVarDbl.h"

using namespace var_dbl;

TestResult test_exp(double value, double uncertainty, std::ostringstream& oss) {
    TestResult res;
    const VarDbl input(value, uncertainty);
    oss << "exp(" << input << ")=";
    res.val = std::exp(value);
    res.expn = VarDbl(std::pow(uncertainty, 2)/2 + std::pow(uncertainty, 4)/8 
                        + std::pow(uncertainty, 6)/48 + std::pow(uncertainty, 8)/384,
                      std::pow(uncertainty, 2) + std::pow(uncertainty, 4)*3/2 
                        + std::pow(uncertainty, 6)*7/6 + std::pow(uncertainty, 8)*5/8, 
                      true);
    std::ostringstream os;
    os << "./Output/exp_" << value << "_" << uncertainty << ".txt";
    res.var = input.exp(os.str().c_str());
    res.err = res.var / res.val - 1;
    return res;
}

VarDbl validate_exp(double value, double uncertainty, double precVal = 1e-3, double precVar = 5e-3)
{
    return validate_func(value, uncertainty, test_exp, precVal, precVar, "");
}

VarDbl validate_exp(double value, double uncertainty, std::string exception)
{
    return validate_func(value, uncertainty, test_exp, 0, 0, exception);
}


int main() 
{
    int i = 19500;
    for (; i < 20000; ++i) {
        try {
            VarDbl(0, i/1000.).exp();
        } catch (const NotMonotonicException& ex) {
            break;
        }
    }
    test::assertEqual(i, 19865);

    std::vector<double> sX{0, 1, -1, 2, -2, 5, -5, 10, -10, 20, -20, 50, -50, 100, -100};
    for (int value: sX) {
        validate_exp(value, 0);
        validate_exp(value, 0.1);
        validate_exp(value, 1, 3e-4, 4e-1);
        validate_exp(value, 19.865, "NotMonotonicException");
        VarDbl res = validate_exp(value, 19.864, 2e+36, 6e+78);
        test::assertAlmostEqual(res.uncertainty() / res.value(), 1680.377, 1e-3);
    }  

    stat_func("Exp", 
            [](double value, double uncertainty){ return VarDbl(value, uncertainty).exp(); },
            [](double value, double noise){ return std::exp(value + noise); },
            sX);
    std::cout << "All exp tests are successful";
}