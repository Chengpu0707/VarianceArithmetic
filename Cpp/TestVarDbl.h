#include "VarDbl.h"
#include "Test.h"
#include "Stat.h"

#ifndef __TestValDbl_h__
#define __TestValDbl_h__
namespace var_dbl 
{

/*
{val}:      the double result of the function   
{var}:      the VarDbl result of the function
{err}:      the difference between var and val
{expn}:     the lower order expansion of err 
*/
struct TestResult {
    VarDbl input, var;
    double val;
    VarDbl err, expn;
};

/*
If {exception} is not empty, expect {func} to throw the corresponding exception.

Otherwise, compare the result {err} and {expn} to the desired {precVal} and {precVar}
*/
VarDbl validate_func(double value, double uncertainty, 
        std::function<TestResult(double, double, std::ostringstream& oss)> func,
        double precVal, double precVar, const std::string exception)
{
    std::ostringstream oss;
    try {
        TestResult res = func(value, uncertainty, oss);
        if (!exception.empty()) {
            oss << res.var;
            test::fail(oss.str());
        }
        oss << res.var << ": " << res.err << " vs " << res.expn;
        test::assertAlmostEqual(res.err.value(), res.expn.value(), precVal, oss.str());
        test::assertAlmostEqual(res.err.variance(), res.expn.variance(), precVar, oss.str());
        return res.var;
    } catch (InitException ex) {
        oss << "InitException";
        test::assertEqual("InitException", exception, oss.str());
        return VarDbl();
    } catch (const TaylorIdException& ex) {
        oss << ex.name();
        test::assertEqual(ex.name(), exception, oss.str());
        return VarDbl();
    }
}

/*
Search the convergence edge for {func}, for x in {sX}, 
    and dx between [sSearch[0]/sSearch[2], sSearch[1]/sSearch[2]]
Output result to file {dumpPath}

If {exception} is not "", {exception} must be thrown for each x in {sX}.
*/
void search_edge(const char* const dumpPath, std::function<VarDbl(double, double)> func,
        std::vector<double> sX, std::array<int, 3> sSearch, bool allowNoEdge,
        const std::string exception = "") 
{
    std::ostringstream oss;
    oss << "x in [";
    for (double x: sX)
        oss << x << ", ";
    oss << "]";
    test::assertLess(0, sX.size(), oss.str());
    oss.str("");
    oss << "search in [";
    for (int i: sSearch)
        oss << i << ",";
    oss << "]";
    test::assertLess(sSearch[0] + 1, sSearch[1], oss.str());
    test::assertLess(0, sSearch[2], oss.str());

    std::ofstream ofs(dumpPath);
    oss.str("");
    oss << "Inalid path " << dumpPath;
    test::assertTrue(ofs.is_open(), oss.str());
    ofs << "X\tEdge\tBias\tValue\tUncertainty\tException\n";

    for (double x: sX) {
        double edge;
        VarDbl res, bias;
        std::string except;
        int i = sSearch[0], j = sSearch[1]; 
        while (i + 1 < j) {
            const int k = (i + j)/2;
            const double dx = ((double) k) / sSearch[2];
            try {
                res = func(x, dx);
                edge = dx;
                i = k;
            } catch (const TaylorIdException& ex) {
                except = ex.name();
                j = k;
            }
        }
        oss.str("");
        oss << "x = " << x << ", i=" << i << ", j=" << i << ", edge=" << edge;
        test::assertNotEqual(i, sSearch[0], oss.str());
        if (allowNoEdge && (j == sSearch[1]))
            continue;
        test::assertNotEqual(j, sSearch[1], oss.str());
        bias = func(x, 0);
        ofs << x << "\t" << edge << "\t" << res.value() - bias.value() << "\t" << res.value() << "\t" << res.uncertainty() << "\t" << except << "\n";
        ofs.flush();
    }
}

/*
Calculate statistics for {test} using {samples} for dx from 1e-1 to 1e-18, for x in {sX}.
    {varFunc} to calculate the expected stat given a value and a deviation
    {dblFunc} to calculate the actual value given a value and a noise
    {samples} is the sample count to construct a stat on the dblFunc using either Gaussian or Uniform sampling
*/
void stat_func(std::string test, 
        std::function<VarDbl(double, double)> varFunc, std::function<double(double, double)> dblFunc,
        std::vector<double> sX, std::vector<double> sDx = std::vector<double>(),
        unsigned samples=10000) 
{
    std::ostringstream oss;
    oss << "x in [";
    for (double x: sX)
        oss << x << ", ";
    oss << "]";
    test::assertLess(0, sX.size(), oss.str());

    if (sDx.empty()) {
        sDx.reserve(18);
        oss.str("");
        oss << "dx=[";
        for (int exp = -17; exp < 0; ++exp) {
            const double dx = std::pow(10, exp);
            sDx.push_back(dx);
            oss << dx << ",";
        }
        oss << "]";
        sDx.push_back(0.2);
        test::assertTrue(std::is_sorted(sDx.begin(), sDx.end()), oss.str());
    }

    Stat<double> stat;
    Histogram histo;

    oss.str("");
    oss << "./Output/" << test << "Stat.txt";
    std::string dumpPath = oss.str();
    std::ofstream ofs(dumpPath);
    oss.str("");
    oss << "Inalid path " << dumpPath;
    test::assertTrue(ofs.is_open(), oss.str());
    ofs << "NoiseType\tNoise\tX\t" << test 
        << "\tError Deviation\tError Minimum\tError Maximum\tValue Deviation\tUncertainty\tMean\tBias"
        << histo.header() << "\n";

    for (double x: sX) {
        for (double dx: sDx) {
            Random rand(0, dx);
            const double dbl = dblFunc(x, 0);
            if (!std::isfinite(dbl)) {
                std::cout << "Ignore " << test << " for x=" << x << " dx=" << dx
                        << " because the result value is infinitive." << std::endl; 
                continue;
            }
            VarDbl v;
            try {
                v = varFunc(x, dx);
                if (v.variance() == 0) {
                    std::cout << "Ignore " << test << " for x=" << x << " dx=" << dx 
                            << " because the result uncertainty is 0." << std::endl; 
                    continue;
                }
            } catch (const std::exception& ex) {
                std::cout << "Ignore " << test << " for x=" << x << " dx=" << dx 
                        << " because " << ex.what() << std::endl; 
                continue;
            }
            const double bias = v.value() - dbl;
            const double unc = v.uncertainty();
            for (bool gauss: {true, false}) {
                stat.clear();
                histo.clear();
                bool finite = true;
                for (size_t i = 0; i < samples; ++i) {
                    const double noise = gauss? rand.gauss() : rand.white();
                    const double err = dblFunc(x, noise) - dbl;
                    if (!std::isfinite(err)) {
                        finite = false;
                        std::cout << "Ignore " << test << " for x=" << x << " dx=" << dx << " gauss=" << gauss << " noise=" << noise 
                                << " because the result contains infinitive after sampling=" << i << "/" << samples << std::endl; 
                        break;
                    }
                    stat.add(err);
                    histo.add(err / unc);
                }
                if (!finite)
                    continue;
                ofs << (gauss? "Gaussian" : "Uniform") 
                    << "\t" << dx << "\t" << x << "\t" << dbl
                    << "\t" << histo.std() << "\t" << histo.min() << "\t" << histo.max()
                    << "\t" << stat.std() << "\t" << unc 
                    << "\t" << stat.mean() << "\t" << bias
                    << histo.formatted() << "\n";
                ofs.flush();
            }
        }
    }

}

}   // namespace var_dbl 
#endif //__TestValDbl_h__