#include "TestTaylor.h"

using namespace var_dbl;

TestResult test_sin(double value, double uncertainty, std::ostringstream& oss) {
    TestResult res;
    res.input = VarDbl(value, uncertainty);
    oss << "sin(" << res.input << ")=";
    res.val = std::sin(value);
    const double cos2 = std::cos(value) * std::cos(value);
    res.expn = VarDbl((std::pow(uncertainty, 2)/2 + std::pow(uncertainty, 4)/8 
                        + std::pow(uncertainty, 6)/48 + std::pow(uncertainty, 8)/384) * res.val,
                       std::sqrt(std::pow(uncertainty, 2) * cos2 - std::pow(uncertainty, 4)*(3./2 * cos2 - 1./2) 
                        + std::pow(uncertainty, 6)*(7./6 * cos2 - 1./2) - std::pow(uncertainty, 8)*(5./8 * cos2 - 7./24)));
    std::ostringstream os;
    os << "./Output/sin_" << value << "_" << uncertainty << ".txt";
    res.var = Taylor::sin(res.input, os.str().c_str());
    res.err = res.var - res.val;
    return res;
}


VarDbl validate_sin(double value, double uncertainty, double precVal = 1e-3, double precVar = 5e-3)
{
    return validate_func(value, uncertainty, test_sin, precVal, precVar, "");
}

VarDbl validate_sin(double value, double uncertainty, std::string exception)
{
    return validate_func(value, uncertainty, test_sin, 0, 0, exception);
}


int main() 
{
    validate_sin(0, 0);
    validate_sin(0, 0.1);
    validate_sin(0, 0.5);
    validate_sin(0, 0.75, 1e-7, 7e-3);
    validate_sin(0, 1, 1e-7, 1e-1);
    validate_sin(0, 1.001, "NotPositiveException");

    validate_sin(0.5 *std::numbers::pi, 0.1, 1e-2);
    validate_sin(0.5 *std::numbers::pi, 1, 1.1, 0.3);

    validate_sin(0.25 *std::numbers::pi, 0.1, 1e-2);
    validate_sin(0.25 *std::numbers::pi, 1, 0.8);

    std::vector<double> sX;
    const int DIVIDS = 64;
    sX.reserve(DIVIDS * 2 + 1);
    for (int j = -DIVIDS; j <= DIVIDS; ++j)
        sX.push_back(std::numbers::pi * j/ DIVIDS);
    search_edge("./Output/SinEdge.txt", 
        [](double value, double uncertainty) { return Taylor::sin(VarDbl(value, uncertainty)); },
        sX, {3000, 50000, 10000}, false, "NotPositiveException");

    stat_func("Sin", 
            [](double value, double uncertainty){ return Taylor::sin(VarDbl(value, uncertainty)); },
            [](double value, double noise){ return std::sin(value + noise); },
            sX);
    std::cout << "All sin tests are successful";
}