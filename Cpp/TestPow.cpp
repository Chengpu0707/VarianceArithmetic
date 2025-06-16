#include "TestVarDbl.h"

using namespace var_dbl;

TestResult test_pow(double value, double uncertainty, std::ostringstream& oss) {
    TestResult res;
    res.input = VarDbl(1, uncertainty);
    oss << "pow(" << res.input << ")^" << value << "=";
    res.val = 1;
    if (0 == value) {
        res.var = VarDbl(1);
        res.err = VarDbl(0, 0); 
        res.expn = VarDbl(0, 0); 
        return res;
    }
    res.expn = VarDbl(std::pow(uncertainty, 2) * value*(value-1)/2 
	                    + std::pow(uncertainty, 4) * value*(value-1)*(value-2)*(value-3)/24
	                    + std::pow(uncertainty, 6) * value*(value-1)*(value-2)*(value-3)*(value-4)*(value-5)/720,
                      std::pow(uncertainty, 2) * value*value
	                    + std::pow(uncertainty, 4) * value*value*(value-1)*(value-5./3)*3./2
	                    + std::pow(uncertainty, 6) * value*value*(value-1)*(value-2)*(value-2)*(value-16./7)*7./6, 
                      true);
    std::ostringstream os;
    os << "./Output/pow_1_" << uncertainty << "_" << value << ".txt";
    res.var = res.input.pow(value, os.str().c_str());
    res.err = res.var - res.val;
    return res;
}

VarDbl validate_pow(double value, double uncertainty, double precVal = 1e-3, double precVar = 5e-3)
{
    return validate_func(value, uncertainty, test_pow, precVal, precVar, "");
}

VarDbl validate_pow(double value, double uncertainty, std::string exception)
{
    return validate_func(value, uncertainty, test_pow, 0, 0, exception);
}


int main() 
{
    VarDbl(1,0.2).pow(-1, "./Output/pow_1_0.2_-2.txt");

    validate_pow(0, 0);
    validate_pow(0, 0.1);
    validate_pow(0, 0.2);
    validate_pow(0, 1);
    validate_pow(0, 10);

    validate_pow(1, 0);
    validate_pow(1, 0.1);
    validate_pow(1, 0.2);
    validate_pow(1, 1);
    validate_pow(1, 10);

    validate_pow(2, 0);
    validate_pow(2, 0.1);
    validate_pow(2, 0.2);
    validate_pow(2, 1);
    validate_pow(2, 10, 2e-3, 4);

    validate_pow(0.5, 0);
    validate_pow(0.5, 0.1);
    validate_pow(0.5, 0.2);
    validate_pow(0.5, 0.21, "NotMonotonicException");

    validate_pow(-1, 0);
    validate_pow(-1, 0.1);
    validate_pow(-1, 0.2, 5e-3, 8e-3);
    validate_pow(-1, 0.21, "NotMonotonicException");

    validate_pow(-2, 0);
    validate_pow(-2, 0.1, 2e-3);
    validate_pow(-2, 0.21, "NotMonotonicException");

    std::vector<double> sX;
    const int DIVIDS = 30;
    sX.reserve(DIVIDS * 2 + 3);
    for (int j = -DIVIDS; j <= DIVIDS; ++j)
        sX.push_back(j /10.);
    search_edge("./Output/PowEdge.txt", 
        [](double value, double uncertainty) { return VarDbl(1, uncertainty).pow(value); },
        sX, {19000, 21000, 100000}, true, "NotMonotonicException");

    sX.push_back(1e-3);
    sX.push_back(-1e-3);
    stat_func("Pow", 
            [](double value, double uncertainty){ return VarDbl(1, uncertainty).pow(value); },
            [](double value, double noise){ return std::pow(1 + noise, value); },
            sX);
    std::cout << "All pow tests are successful";
}