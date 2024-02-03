#include "VarDbl.h"
#include "Test.h"

#include <numbers>

using namespace var_dbl;

void assertEquals(VarDbl var, double value, double uncertainty) 
{
    test::assertEquals(var.value(), value);
    test::assertEquals(var.uncertainty(), uncertainty);
}

constexpr static const double EXP_VALUE_DELTA = 5e-7;
constexpr static const double EXP_VARIANCE_DELTA = 1e-6;

constexpr static const double SIN_VALUE_DELTA = 5e-7;
constexpr static const double SIN_VARIANCE_DELTA = 1e-6;

constexpr static const double LOG_VALUE_DELTA = 2e-3;
constexpr static const double LOG_VARIANCE_DELTA = 5e-3;

constexpr static const double POW_VALUE_DELTA = 2e-3;
constexpr static const double POW_VARIANCE_DELTA = 5e-3;



VarDbl validate_exp(double exp, double uncertainty, 
                    double valueDelta=EXP_VALUE_DELTA, double varianceDelta=EXP_VARIANCE_DELTA)
{
    const VarDbl var(exp, uncertainty);
    const VarDbl res = var.exp();
    const VarDbl prec = res * std::exp(-exp);
    std::ostringstream os;
    os << "exp(" << var << ")=" << res;
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


VarDbl validate_log(double x, double uncertainty, 
                    double valueDelta=LOG_VALUE_DELTA, double varianceDelta=LOG_VARIANCE_DELTA)
{
    const VarDbl var(x, uncertainty);
    const VarDbl res = var.log();
    const VarDbl prec = res - std::log(x);
    std::ostringstream os;
    os << "log(" << var << ")=" << res;
    const double precIn = abs(uncertainty/x);
    test::assertEquals(prec.value(),
            - std::pow(precIn, 2)/2 - std::pow(precIn, 4)/4 \
            - std::pow(precIn, 6)/6 - std::pow(precIn, 8)/8,
            valueDelta, os.str());
    test::assertEquals(prec.variance(),
            std::pow(precIn, 2) + std::pow(precIn, 4)*9/8 +\
            std::pow(precIn, 6)*119/24 + std::pow(precIn, 8)*991/32,
            varianceDelta, os.str());
    return res;
}


void test_taylor_log(std::initializer_list<double> sValue, double uncertainty, 
                     double valueDelta=LOG_VALUE_DELTA, double varianceDelta=LOG_VARIANCE_DELTA)
{
    try {
        validate_log(*sValue.begin(), uncertainty);
    } catch (LossUncertaintyException ex) {
    }
    try {
        validate_log(*(sValue.begin() + 1), uncertainty);
    } catch (test::AssertException ex) {
    }
    validate_log(*(sValue.begin() + 2), uncertainty, valueDelta, varianceDelta);

}

void test_taylor_log()
{   
    test_taylor_log({1./4, 1./2, 1.}, 0.2);
    test_taylor_log({1./8, 1./4, 1./2}, 0.1);
    test_taylor_log({1./64, 1./32, 1./16}, 0.01);
}

VarDbl validate_sin(double x, double uncertainty, 
                    double valueDelta=SIN_VALUE_DELTA, double varianceDelta=SIN_VARIANCE_DELTA)
{
    const VarDbl var(x, uncertainty);
    const VarDbl res = var.sin();
    const VarDbl prec = res - std::sin(x);
    std::ostringstream os;
    os << "sin(" << var << ")=" << res;
    test::assertEquals(prec.value(),
            std::sin(x) * (-std::pow(uncertainty, 2)/2 + std::pow(uncertainty, 4)/8 +\
                           -std::pow(uncertainty, 6)/48 + std::pow(uncertainty, 8)/384),
            valueDelta, os.str());
    const double cos2 = std::cos(x) * std::cos(x);
    test::assertEquals(prec.variance(),
            std::pow(uncertainty, 2) * cos2 - std::pow(uncertainty, 4)*(3./2 * cos2 - 1./2) +\
                        std::pow(uncertainty, 6)*(7./6 * cos2 - 1./2),
            varianceDelta, os.str());
    return res;
}

void test_taylor_sin()
{   
    const double PI = 3.14159265358979323846;
    validate_sin(0, 0);
    validate_sin(0, 1e-3);
    validate_sin(0, 1e-2);
    validate_sin(0, 1e-1);

    validate_sin(PI/2, 0);
    validate_sin(PI/2, 1e-3);
    validate_sin(PI/2, 1e-2);
    validate_sin(PI/2, 1e-1);

    validate_sin(-PI/2, 0);
    validate_sin(-PI/2, 1e-3);
    validate_sin(-PI/2, 1e-2);
    validate_sin(-PI/2, 1e-1);
}


VarDbl validate_pow(double exp, double uncertainty, 
                    double valueDelta=POW_VALUE_DELTA, double varianceDelta=POW_VARIANCE_DELTA)
{
    const VarDbl var(1, uncertainty);
    const VarDbl res = var.pow(exp);
    const VarDbl prec = res - VarDbl(1, 0);
    std::ostringstream os;
    os << "pow(" << var << ", " << exp << ")=" << res;
    test::assertEquals(prec.value(),
	    std::pow(uncertainty, 2) * exp*(exp-1)/2 + \
	    std::pow(uncertainty, 4) * exp*(exp-1)*(exp-2)*(exp-3)/24 +\
	    std::pow(uncertainty, 6) * exp*(exp-1)*(exp-2)*(exp-3)*(exp-4)*(exp-5)/720,
            valueDelta, os.str());
    test::assertEquals(prec.variance(),
	    std::pow(uncertainty, 2) * exp*exp + \
	    std::pow(uncertainty, 4) * exp*exp*(exp-1)*(exp-5./3)*3./2 +\
	    std::pow(uncertainty, 6) * exp*exp*(exp-1)*(exp-2)*(exp-2)*(exp-16./7)*7./6,
            varianceDelta, os.str());
    return res;
}


void test_taylor_pow()
{
    validate_pow(0, 0); 
    validate_pow(0, 0.2); 
    validate_pow(1e-6, 0);    
    validate_pow(1e-3, 0);   
    validate_pow(0.1, 0);   

    validate_pow(1, 0.2);  
    validate_pow(2, 0.2);    
    validate_pow(0.5, 0.2);    
    validate_pow(0.5, 0.205);    
    validate_pow(0.5, 0.23,    7e-3, 2e-2);    
    validate_pow(-1, 0.2,      5e-3, 7e-3); 
    validate_pow(-1, 0.195,    5e-3, 7e-3);  
    validate_pow(-1.5, 0.205,  2e-2, 0.6); 
    validate_pow(-1.5, 0.2,    2e-2, 8e-2); 
    validate_pow(-1.5, 0.195,  2e-2, 3e-2);  
    validate_pow(-2, 0.2,      3e-2, 1.7);  
    validate_pow(-2, 0.195,    3e-2, 3e-1);  
}


int main() 
{
    test_taylor_exp();
    test_taylor_log();
    test_taylor_sin();
    test_taylor_pow();
    std::cout << "All Taylor expansion tests are successful";
    return 0;
}