#include "TestVarDbl.h"

using namespace var_dbl;

TestResult test_log(double value, double uncertainty, std::ostringstream& oss) {
    TestResult res;
    res.input = VarDbl(value, uncertainty);
    oss << "log(" << res.input << ")=";
    res.val = std::log(value);
    const double precIn = uncertainty / value;
    res.expn = VarDbl(- std::pow(precIn, 2)/2 - std::pow(precIn, 4)/4 
                        - std::pow(precIn, 6)/6 - std::pow(precIn, 8)/8,
                      std::pow(precIn, 2) + std::pow(precIn, 4)*9/8 
                        + std::pow(precIn, 6)*119/24 + std::pow(precIn, 8)*991/32, 
                      true);
    std::ostringstream os;
    os << "./Output/log_" << value << "_" << uncertainty << ".txt";
    res.var = res.input.log(os.str().c_str());
    res.err = res.var - res.val;
    return res;
}

VarDbl validate_log(double value, double uncertainty, double precVal = 1e-3, double precVar = 5e-3)
{
    return validate_func(value, uncertainty, test_log, precVal, precVar, "");
}

VarDbl validate_log(double value, double uncertainty, std::string exception)
{
    return validate_func(value, uncertainty, test_log, 0, 0, exception);
}


int main() 
{
    int i = 20080;
    for (; i < 20100; ++i) {
        try {
            VarDbl(1, i/100000.).log();
        } catch (std::runtime_error ex) {
            break;
        }
    }
    test::assertEqual(i, 20087);
    std::vector<double> sX{1., 2., 1./2, 4., 1./4, 8., 1./8, 16., 1./16, 64., 1./64};
    for (double value: sX) {
        validate_log(value, 0);
        validate_log(value, 0.1 * value);
        validate_log(value, 0.15 * value);
        validate_log(value, 0.20087 * value, "NotMonotonicException");
        const VarDbl res = validate_log(value, 0.20086 * value, 2e-3);
        test::assertAlmostEqual(res.uncertainty(), 0.2130506, 1e-6);
    } 

    stat_func("Log", 
            [](double value, double uncertainty){ return VarDbl(value, uncertainty).log(); },
            [](double value, double noise){ return std::log(value + noise); },
            sX);
    std::cout << "All log tests are successful";
}