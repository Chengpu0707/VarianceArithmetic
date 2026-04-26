#include "TestTaylor.h"

#if __cplusplus >= 202002L
#include <numbers>
#define _TEST_SIN_PI std::numbers::pi
#else
static const double _test_sin_pi = 3.14159265358979323846;
#define _TEST_SIN_PI _test_sin_pi
#endif

using namespace var_dbl;

#if __cplusplus < 201103L
struct _SinEdgeFunctor {
    VarDbl operator()(double value, double uncertainty) const {
        return Taylor::sin(VarDbl(value, uncertainty));
    }
};
struct _SinDblFunctor {
    double operator()(double value, double noise) const {
        return std::sin(value + noise);
    }
};
#endif

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

    validate_sin(0.5 *_TEST_SIN_PI, 0.1, 1e-2);
    validate_sin(0.5 *_TEST_SIN_PI, 1, 1.1, 0.3);

    validate_sin(0.25 *_TEST_SIN_PI, 0.1, 1e-2);
    validate_sin(0.25 *_TEST_SIN_PI, 1, 0.8);

    std::vector<double> sX;
    const int DIVIDS = 64;
    sX.reserve(DIVIDS * 2 + 1);
    for (int j = -DIVIDS; j <= DIVIDS; ++j)
        sX.push_back(_TEST_SIN_PI * j/ DIVIDS);
#if __cplusplus >= 201103L
    search_edge("./Output/SinEdge.txt",
        [](double value, double uncertainty) { return Taylor::sin(VarDbl(value, uncertainty)); },
        sX, {3000, 50000, 10000}, false, "NotPositiveException");
#else
    { int _sSearch[] = {3000, 50000, 10000};
      search_edge("./Output/SinEdge.txt", _SinEdgeFunctor(), sX, _sSearch, false, "NotPositiveException"); }
#endif

#if __cplusplus >= 201103L
    stat_func("Sin",
            [](double value, double uncertainty){ return Taylor::sin(VarDbl(value, uncertainty)); },
            [](double value, double noise){ return std::sin(value + noise); },
            sX);
#else
    stat_func("Sin", _SinEdgeFunctor(), _SinDblFunctor(), sX);
#endif
    std::cout << "All sin tests are successful";
}