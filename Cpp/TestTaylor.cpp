#include "VarDbl.h"
#include "Test.h"
#include "ulp.h"

using namespace var_dbl;

void assertEquals(VarDbl var, double value, double uncertainty) 
{
    test::assertEquals(var.value(), value);
    test::assertEquals(var.uncertainty(), uncertainty);
}


constexpr static const double EXP_VALUE_DELTA = 5e-7;
constexpr static const double EXP_VARIANCE_DELTA = 1e-6;

VarDbl validate_exp(double exp, double uncertainty, 
                    double valueDelta=EXP_VALUE_DELTA, double varianceDelta=EXP_VARIANCE_DELTA)
{
    const VarDbl var(exp, uncertainty);
    const VarDbl res = var.exp();
    const VarDbl prec = res * std::exp(-exp);
    std::ostringstream os;
    os << var << ".exp()=" << res;
    test::assertEquals(prec.value(),
            1 + std::pow(uncertainty, 2)/2 + std::pow(uncertainty, 4)/8 +\
                std::pow(uncertainty, 6)/48 + std::pow(uncertainty, 8)/384,
            valueDelta, os.str());
    test::assertEquals(prec.variance(),
            std::pow(uncertainty, 2) + std::pow(uncertainty, 4)*3/2 +\
                std::pow(uncertainty, 6)*7/6 + std::pow(uncertainty, 8)*5/8,
            varianceDelta, os.str());
    return res;
}


void test_taylor_exp(std::initializer_list<double> sValue, double uncertainty, 
                     double valueDelta=EXP_VALUE_DELTA, double varianceDelta=EXP_VARIANCE_DELTA)
{
    try {
        validate_exp(*sValue.begin(), uncertainty);
    } catch (UncertaintyException ex) {
    }
    try {
        validate_exp(*(sValue.begin() + 1), uncertainty);
    } catch (test::AssertException ex) {
    }
    validate_exp(*(sValue.begin() + 2), uncertainty, valueDelta, varianceDelta);

}

void test_taylor_exp()
{   
    test_taylor_exp({-392., -390., -200.}, 0.1);
    test_taylor_exp({196., 194., 100.}, 0.1);

    test_taylor_exp({-392., -390., -200.}, 0.2, 1e-3);
    test_taylor_exp({196., 194., 100.}, 0.2, 1e-3);

    test_taylor_exp({-392., -390., -200.}, 1, 1e-3, 0.5);
    test_taylor_exp({196., 194., 100.}, 1, 1e-3, 0.5);
}




int main() 
{
    test_taylor_exp();

    std::cout << "All Taylor expansion tests are successful";
    return 0;
}