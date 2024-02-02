
#include "VarDbl.h"
#include "Stat.h"
#include "Test.h"

#include <functional>
#include <fstream>
#include <random>
#include <sstream>
#include <vector>

using namespace var_dbl;

void testGuassian()
{
    const size_t SAMPLES = 1000;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution normal{0., 1.};
    double sNoise[SAMPLES];
    for (int i = 0; i < SAMPLES; ++i)
        sNoise[i] = normal(gen) * 0.2 + 1;
    Stat stat = calcStat(sNoise, sNoise + SAMPLES);
    test::assertEquals(stat.mean, 1, 1e-2);
    test::assertEquals(stat.stddev, 0.2, 1e-1);
}

void testUniform() 
{
    const size_t SAMPLES = 1000;
    double sUniform[SAMPLES];
    int half = SAMPLES / 2;
    for (int i = -half; i < half; ++i)
        sUniform[i + half] = std::sqrt(3) *i /half *0.2 + 1;
    Stat stat = calcStat(sUniform, sUniform + SAMPLES);
    test::assertEquals(stat.mean, 1, 1e-2);
    test::assertEquals(stat.stddev, 0.2, 1e-1);
}


bool dump_test(const std::string test, const std::vector<double>& sDev, const std::vector<double>& sX)
{
    const float HIST_RANGE = 3;
    const unsigned HIST_DIVIDES = 5;
    const unsigned SAMPLES = 10000;

    std::ostringstream os;
    os << "./Output/" << test << "Var.txt";
    std::ofstream of(os.str());
    of << "NoiseType\tNoise\tX\t" << test << "\tError Deviation\tError Minimum\tError Maximum"
        << "\tValue Deviation\tUncertainty\tMean\tBias";
    auto sHistCenter = calcHistoCenters(HIST_RANGE, HIST_DIVIDES);
    for (auto center: sHistCenter)
        of << "\t" << center;
    of << "\n";

    double sNormal[SAMPLES], sUniform[SAMPLES];
    double sValueError[SAMPLES], sNormError[SAMPLES];
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution normal{0., 1.};
    for (int i = 0; i < SAMPLES; ++i)
        sNormal[i] = normal(gen);
    int half = SAMPLES / 2;
    for (int i = -half; i < half; ++i)
        sUniform[i + half] = std::sqrt(3) *i /half;
    Histo histo;

    for (auto dev: sDev) {
        for (auto x: sX) {
            const VarDbl var(x, dev);
            double base;
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
                const double* sNoise = (iNoise == 0)? sNormal : sUniform;
                bool notFinite = false;
                for (size_t i = 0; i < SAMPLES; ++i) {
                    if (test == "exp")
                        sValueError[i] = std::exp(x + sNoise[i] * dev) - base;
                    else if (test == "log")
                        sValueError[i] = std::log(x + sNoise[i] * dev) - base;
                    else if (test == "sin")
                        sValueError[i] = std::sin(x + sNoise[i] * dev) - base;
                    else
                        return false;
                    if (!std::isfinite(sValueError[i])) {
                        std::cout << "Ignore " << test << "(" << var << ")=" << res 
                            <<": result is not finite for iNoise="<< iNoise << "\n";
                        notFinite = true;
                        break;
                    }
                    sNormError[i] = sValueError[i] / res.uncertainty();
                }
                if (notFinite)
                    continue;
                Stat stat = calcStat(sValueError, sValueError + SAMPLES);
                if (!calcHisto(histo, sNormError, sNormError + SAMPLES, HIST_RANGE, HIST_DIVIDES)) {
                    std::cout << "Ignore " << test << " of " << var << "=" << res 
                        <<": histogram fails for iNoise="<< iNoise << "\n";
                    continue;
                }
                of << ((iNoise == 0)? "Gaussian" : "Uniform") 
                    << "\t" << dev << "\t" << x << "\t" << base
                    << "\t" << histo.stat.stddev << "\t" << histo.stat.min << "\t" << histo.stat.max
                    << "\t" << stat.stddev << "\t" << res.uncertainty() 
                    << "\t" << stat.mean << "\t" << res.value() - base;
                for (auto prob: histo.sHisto) 
                    of << "\t" << prob;
                of << "\n";
            }
        }
    }
    std::cout << "finish write " << os.str();
    return true;
}            




int main()
{
    testGuassian();
    testUniform();

    const std::vector<double> DEVS{
        0.2, 1e-2, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8,
        1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17
    };
    std::vector<double> sDev(DEVS);
    sDev[0] = 1;

//    test::assertTrue( dump_test("exp", sDev, {-100, -50, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 50, 100}) );

//    test::assertTrue( dump_test("log", DEVS, {1./32, 1./20, 1./16, 0.1, 1./8, 1./4, 1./2, 1, 2, 4, 8, 16, 32}) );
    const double PI = 3.14159265358979323846;
    std::vector<double> sSinX;
    for (int i = -16; i <= 16; ++i) 
        sSinX.push_back(PI*i/16);
    test::assertTrue( dump_test("sin", sDev, sSinX) );
}