#include "VarDbl.h"
#include "Test.h"
#include "Stat.h"

#include <chrono>
#include <fstream>
#include <numbers>

using namespace var_dbl;

void assertEquals(VarDbl var, double value, double uncertainty) 
{
    test::assertAlmostEqual(var.value(), value);
    test::assertAlmostEqual(var.uncertainty(), uncertainty);
}

constexpr static const double VALUE_DELTA = 1e-3;
constexpr static const double VARIANCE_DELTA = 5e-3;

constexpr const char* EDGE_HEADER = "Edge Value\tEdge Uncertainty\tValue\tUncertainty\tExceptionType\tOrder\n";


VarDbl validate_exp(double exp, double uncertainty)
{
    const VarDbl var(exp, uncertainty);
    const VarDbl res = var.exp();
    const VarDbl prec = res * std::exp(-exp);
    std::ostringstream os;
    os << "exp(" << var << ")=" << res;
    test::assertAlmostEqual(prec.value(),
            1 + std::pow(uncertainty, 2)/2 + std::pow(uncertainty, 4)/8 +\
                std::pow(uncertainty, 6)/48 + std::pow(uncertainty, 8)/384,
            VALUE_DELTA, os.str());
    test::assertAlmostEqual(prec.variance(),
            std::pow(uncertainty, 2) + std::pow(uncertainty, 4)*3/2 +\
                std::pow(uncertainty, 6)*7/6 + std::pow(uncertainty, 8)*5/8,
            VARIANCE_DELTA, os.str());
    return res;
}


void test_taylor_exp(std::initializer_list<double> sValue, double uncertainty)
{
    try {
        validate_exp(*sValue.begin(), uncertainty);
        test::fail("not detect infinity");
    } catch (InitException ex) {
    }
    validate_exp(*(sValue.begin() + 1), uncertainty);
}

void test_taylor_exp()
{   
    test_taylor_exp({-392., -370.}, 0.1);
    test_taylor_exp({-392., -370.}, 0.2);
    test_taylor_exp({-392., -370.}, 0.5);
    test_taylor_exp({-392., 0.}, 0.66);

    test_taylor_exp({195.671, 195.670}, 0.1);
    test_taylor_exp({195.671, 195.670}, 0.2);
    test_taylor_exp({195.671, 195.670}, 0.5);
    test_taylor_exp({195.671, 190.670}, 0.66);
}


VarDbl validate_log(double x, double uncertainty)
{
    const VarDbl var(x, uncertainty);
    const VarDbl res = var.log();
    const VarDbl prec = res - std::log(x);
    std::ostringstream os;
    os << "log(" << var << ")=" << res;
    const double precIn = abs(uncertainty/x);
    test::assertAlmostEqual(prec.value(),
            - std::pow(precIn, 2)/2 - std::pow(precIn, 4)/4 \
            - std::pow(precIn, 6)/6 - std::pow(precIn, 8)/8,
            VALUE_DELTA, os.str());
    test::assertAlmostEqual(prec.variance(),
            std::pow(precIn, 2) + std::pow(precIn, 4)*9/8 +\
            std::pow(precIn, 6)*119/24 + std::pow(precIn, 8)*991/32,
            VARIANCE_DELTA, os.str());
    return res;
}


void test_taylor_log(std::initializer_list<double> sValue, double uncertainty)
{
    try {
        validate_log(*sValue.begin(), uncertainty);
        test::fail("not detect not monotonic");
    } catch (NotMonotonicException ex) {
    }
    validate_log(*(sValue.begin() + 1), uncertainty);
}

void test_taylor_log()
{   
    test_taylor_log({0.95, 1.}, 0.2);
    test_taylor_log({0.493, 0.5}, 0.1);
    test_taylor_log({0.0493, 0.05}, 0.01);
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
                    double valueDelta=VALUE_DELTA, double varianceDelta=VARIANCE_DELTA)
{
    const VarDbl var(x, uncertainty);
    const VarDbl res = var.sin();
    const VarDbl prec = res - std::sin(x);
    std::ostringstream os;
    os << "sin(" << var << ")=" << res;
    test::assertAlmostEqual(prec.value(),
            std::sin(x) * (-std::pow(uncertainty, 2)/2 + std::pow(uncertainty, 4)/8 +\
                           -std::pow(uncertainty, 6)/48 + std::pow(uncertainty, 8)/384),
            valueDelta, os.str());
    const double cos2 = std::cos(x) * std::cos(x);
    test::assertAlmostEqual(prec.variance(),
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
    validate_sin(0, std::numbers::pi*10/32, VALUE_DELTA, 3e-1);
    try {
        const VarDbl var(0, std::numbers::pi*11/32);
        const VarDbl res = var.sin();
        test::fail("sin(0~std::numbers::pi*11/32)");
    } catch (NotReliableException ex) {
    }
}

void test_taylor_sin_half_pi()
{   
    validate_sin(std::numbers::pi/2, 0);
    validate_sin(std::numbers::pi/2, 1e-3);
    validate_sin(std::numbers::pi/2, 1e-2);
    validate_sin(std::numbers::pi/2, 1e-1);
    validate_sin(std::numbers::pi/2, std::numbers::pi*10/32, 3e-4, 2e-1);
    try {
        const VarDbl var(std::numbers::pi/2, std::numbers::pi*11/32);
        const VarDbl res = var.sin();
        test::fail("sin(std::numbers::pi/2~std::numbers::pi*11/32)");
    } catch (NotReliableException ex) {
    }
}

void test_taylor_sin_quarter_pi()
{   
    validate_sin(std::numbers::pi/4, 0);
    validate_sin(std::numbers::pi/4, 1e-3);
    validate_sin(std::numbers::pi/4, 1e-2);
    validate_sin(std::numbers::pi/4, 1e-1);
    validate_sin(std::numbers::pi/4, 1.41, 5e-3, 3e-1);
    try {
        const VarDbl var(std::numbers::pi/4, 1.42);
        const VarDbl res = var.sin();
        test::fail("sin(std::numbers::pi/4~1.41)");
    } catch (NotReliableException ex) {
    }
}


VarDbl validate_pow(double exp, double uncertainty, 
                    double valueDelta=VALUE_DELTA, double varianceDelta=VARIANCE_DELTA)
{
    const VarDbl var(1, uncertainty);
    const VarDbl res = var.pow(exp);
    const VarDbl prec = res - VarDbl(1, 0);
    std::ostringstream os;
    os << "pow(" << var << ", " << exp << ")=" << res;
    test::assertAlmostEqual(prec.value(),
	    std::pow(uncertainty, 2) * exp*(exp-1)/2 + \
	    std::pow(uncertainty, 4) * exp*(exp-1)*(exp-2)*(exp-3)/24 +\
	    std::pow(uncertainty, 6) * exp*(exp-1)*(exp-2)*(exp-3)*(exp-4)*(exp-5)/720,
            valueDelta, os.str());
    test::assertAlmostEqual(prec.variance(),
	    std::pow(uncertainty, 2) * exp*exp + \
	    std::pow(uncertainty, 4) * exp*exp*(exp-1)*(exp-5./3)*3./2 +\
	    std::pow(uncertainty, 6) * exp*exp*(exp-1)*(exp-2)*(exp-2)*(exp-16./7)*7./6,
            varianceDelta, os.str());
    return res;
}


void dump_power_at_one_edge()
    // 97 seconds
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


void dump_sin_edge()
    // 97 seconds
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::ofstream of("./Output/SinEdge.txt");
    of << EDGE_HEADER;
    for (int i = -16; i <= 16; ++i) {
        const double x = i/16. * std::numbers::pi;
        int order = 0;
        std::string type;
        for (int j = 150; j >= 50; --j) {
            const double dx = j/100.;
            const VarDbl var(x, dx);
            try {
                const VarDbl res = var.sin();
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
    std::cout << "dump_sin_at_one_edge() = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[sec]" << std::endl;
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


bool dump_test(const std::string test, const std::vector<double>& sDev, const std::vector<double>& sX)
{
    const unsigned SAMPLES = 10000;
    Stat<double> stat;
    Histogram histo;

    std::ostringstream os;
    os << "./Output/" << test << "Var.txt";
    std::ofstream of(os.str());
    of << "NoiseType\tNoise\tX\t" << test << "\tError Deviation\tError Minimum\tError Maximum"
        << "\tValue Deviation\tUncertainty\tMean\tBias"
        << histo.header() << "\n";

    for (auto dev: sDev) {
        for (auto x: sX) {
            const VarDbl var((test == "pow")? 1 : x, dev);
            double base;
            Random rand(0, dev);
            VarDbl res;
            try {
                if (test == "exp") {
                    base = std::exp(x);
                    res = var.exp();
                } else if (test == "log") {
                    base = std::log(x);
                    res = var.log();
                } else if (test == "sin") {
                    base = std::sin(x);
                    res = var.sin();
                } else if (test == "pow") {
                    base = 1;
                    res = var.pow(x);
                } else
                    return false;
            } catch (std::runtime_error ex) {
                std::cout << "Ignore " << test << "(" << var << "): " << ex.what() << "\n";
                continue;
            }
            if (res.uncertainty() <= 0) {
                std::cout << "Ignore " << test << "(" << var << "): 0=uncertainty for " << res << "\n";
                continue;
            }
            for (size_t iNoise = 0; iNoise < 2; ++iNoise) {
                stat.clear();
                histo.clear();
                bool notFinite = false;
                for (size_t i = 0; i < SAMPLES; ++i) {
                    const double noise = (iNoise == 0)? rand.gauss() : rand.white();
                    double valErr, nrmErr;
                    if (test == "exp")
                        valErr = std::exp(x + noise) - base;
                    else if (test == "log")
                        valErr = std::log(x + noise) - base;
                    else if (test == "sin")
                        valErr = std::sin(x + noise) - base;
                    else if (test == "pow")
                        valErr = std::pow(1 + noise, x) - base;
                    else
                        return false;
                    if (!std::isfinite(valErr)) {
                        std::cout << "Ignore " << test << "(" << var << ")=" << res 
                            <<": result is not finite for iNoise="<< iNoise << "\n";
                        notFinite = true;
                        break;
                    }
                    nrmErr = valErr / res.uncertainty();
                    stat.add(valErr);
                    histo.add(nrmErr);
                }
                if (notFinite)
                    continue;
                of << ((iNoise == 0)? "Gaussian" : "Uniform") 
                    << "\t" << dev << "\t" << x << "\t" << base
                    << "\t" << histo.std() << "\t" << histo.min() << "\t" << histo.max()
                    << "\t" << stat.std() << "\t" << res.uncertainty() 
                    << "\t" << stat.mean() << "\t" << res.value() - base
                    << histo.formatted() << "\n";
            }
        }
    }
    std::cout << "finish write " << os.str() << "\n";
    return true;
}            




void dump_taylor()
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::vector<double> DEVS;
    for (int i = -17; i <= 0; ++i) 
        DEVS.push_back(std::pow(10, i));
    std::vector<double> sDev(DEVS);
    sDev[sDev.size() -1] = 0.2;

    test::assertTrue( dump_test("exp", DEVS, 
            {0, 1, -1, 2, -2, 5, -5, 10, -10, 20, -20, 50, -50, 100, -100}) );

    test::assertTrue( dump_test("log", sDev, 
            {1, 2, 1./2, 4, 1./4, 8, 1./8, 10, 0.1, 16, 1./16, 32, 1./32}) );

    std::vector<double> sSinX;
    for (int i = -16; i <= 16; ++i) 
        sSinX.push_back(std::numbers::pi*i/16);
    test::assertTrue( dump_test("sin", sDev, sSinX) );

    sDev[0] = 0.2;
    sDev.insert(sDev.begin() + 1, 0.195);
    std::vector<double> sPowX;
    for (int i = -20; i <= 32; i += 2) 
        if (i != 0)
            sPowX.push_back(i/10.);
    sPowX.insert(sPowX.end(), {-1e-1, 1e-1, -1e-2, 1e-2, -1e-3, 1e-3, -1e-6, 1e-6});
    test::assertTrue( dump_test("pow", sDev, sPowX) );
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "dump_taylor() = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[sec]" << std::endl;
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
    // dump_power_at_one_edge();   // 94 sec
    dump_sin_edge();    // 0 sec

    dump_taylor();  // 14 sec

    std::cout << "All Taylor expansion tests are successful";
    return 0;
}