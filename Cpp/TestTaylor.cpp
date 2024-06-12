#include "VarDbl.h"
#include "Test.h"

#include <chrono>
#include <fstream>
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

constexpr static const double PI = 3.141592653589793238463;

constexpr const char* EDGE_HEADER = "Edge Value\tEdge Uncertainty\tValue\tUncertainty\tExceptionType\tOrder\n";


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
    } catch (InitException ex) {
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
    } catch (NotMonotonicException ex) {
    }
    try {
        validate_log(*(sValue.begin() + 1), uncertainty);
    } catch (test::AssertException ex) {
    }
    validate_log(*(sValue.begin() + 2), uncertainty, valueDelta, varianceDelta);

}

void test_taylor_log()
{   
    test_taylor_log({15./16, 255./256, 1.}, 0.2);
    test_taylor_log({15./32, 127./256, 1./2}, 0.1);
    test_taylor_log({0.05 - 1e-5, 0.05 - 1e-6, 0.05}, 0.01);
}

void dump_log_edge()
    // 7 seconds
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::ofstream of("./Output/LogEdge.txt");
    of << EDGE_HEADER;
    for (int i = -15; i <= 16; i += 1) {
        const double x = i / 16. + 1;
        double intpart;
        int order = 0;
        std::string type;
        for (int j = 300; j > 0; j -= 1) {
            const double dx = j / 1000.;
            const VarDbl var(x, dx);
            try {
                const VarDbl res = var.log();
                of << x << "\t" << dx << "\t"  << res.value() << "\t" << res.uncertainty()
                   << "\t" << type << "\t" << order << "\n";
                of.flush();
                break;
            } catch (DivergentException ex) {
                type = "DivergentException";
                order = ex.n;
            } catch (NotMonotonicException ex) {
                type = "NotMonotonicException";
                order = ex.n;
            } catch (NotReliableException ex) {
                type = "NotReliableException";
                order = ex.n;
            }
        }
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "dump_log_at_one_edge() = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[sec]" << std::endl;
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

void test_taylor_sin_pi()
{   
    validate_sin(0, 0);
    validate_sin(0, 1e-3);
    validate_sin(0, 1e-2);
    validate_sin(0, 1e-1);
    validate_sin(0, PI*10/32, SIN_VALUE_DELTA, 3e-1);
    try {
        const VarDbl var(0, PI*11/32);
        const VarDbl res = var.sin();
        test::fail("sin(0~PI*11/32)");
    } catch (NotReliableException ex) {
    }
}

void test_taylor_sin_half_pi()
{   
    validate_sin(PI/2, 0);
    validate_sin(PI/2, 1e-3);
    validate_sin(PI/2, 1e-2);
    validate_sin(PI/2, 1e-1);
    validate_sin(PI/2, PI*10/32, 3e-4, 2e-1);
    try {
        const VarDbl var(PI/2, PI*11/32);
        const VarDbl res = var.sin();
        test::fail("sin(PI/2~PI*11/32)");
    } catch (NotReliableException ex) {
    }
}

void test_taylor_sin_quarter_pi()
{   
    validate_sin(PI/4, 0);
    validate_sin(PI/4, 1e-3);
    validate_sin(PI/4, 1e-2);
    validate_sin(PI/4, 1e-1);
    validate_sin(PI/4, 1.41, 5e-3, 3e-1);
    try {
        const VarDbl var(PI/4, 1.42);
        const VarDbl res = var.sin();
        test::fail("sin(PI/4~1.41)");
    } catch (NotReliableException ex) {
    }
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


void dump_power_at_one_edge()
    // 128 seconds
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::ofstream of("./Output/PowerAtOneEdge.txt");
    of << EDGE_HEADER;
    for (int i = -600; i <= 600; i += 2) {
        if ((i % 100) == 0.)
            continue;
        const double exp = i / 100.;
        int order = 0;
        std::string type;
        for (int j = 300; j > 0; j -= 1) {
            const double dx = j / 1000.;
            const VarDbl var(1., dx);
            try {
                const VarDbl res = var.pow(exp);
                of << exp << "\t" << dx << "\t"  << res.value() << "\t" << res.uncertainty()
                   << "\t" << type << "\t" << order << "\n";
                of.flush();
                break;
            } catch (DivergentException ex) {
                type = "DivergentException";
                order = ex.n;
            } catch (NotMonotonicException ex) {
                type = "NotMonotonicException";
                order = ex.n;
            } catch (NotReliableException ex) {
                type = "NotReliableException";
                order = ex.n;
            }
        }
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "dump_power_at_one_edge() = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[sec]" << std::endl;
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
    validate_pow(0.5, 0.203);    
    validate_pow(0.5, 0.2);    
    validate_pow(-1, 0.199,   5e-3, 7e-3); 
    validate_pow(-1.5, 0.197, 2e-2, 0.6); 
    validate_pow(-2, 0.194,   3e-2, 3e-1);  
}


int main() 
{
    test_taylor_exp();
    test_taylor_log();
    test_taylor_sin_pi();
    test_taylor_sin_half_pi();
    test_taylor_sin_quarter_pi();
    test_taylor_pow();

    dump_log_edge();    // 4 sec
    dump_power_at_one_edge();   // 94 sec

    std::cout << "All Taylor expansion tests are successful";
    return 0;
}