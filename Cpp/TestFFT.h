#include <map>
#if __cplusplus >= 201103L
#include <array>
#include <string>
#endif

#include "FFT.h"

#ifndef __Test_FFT_h__
#define __Test_FFT_h__
namespace var_dbl
{

/*
Signals for FFT calculation as sWave and sFreq for different order, frequence and sin source (index vs library).
The Aggr signal is the average result of All Sin and Cos, so that it is not a signal per se.
*/
class FFT_Signal : public FFT {
public:
    enum SignalType {
        Sin,
        Cos,
        Linear,
        Aggr,  // aggregated Sin and Cos for all frequencies
    };
#if __cplusplus >= 202002L
    constexpr static const std::array<std::string, 4> sSignalType{ "Sin", "Cos", "Linear", "Aggr" };
    constexpr static const std::string signalTypeName(SignalType signalType) { return sSignalType[static_cast<size_t>(signalType)]; }
#if __cplusplus >= 202302L
    constexpr
#endif
    static SignalType toSignalType(const std::string & name) {
        for (size_t i = 0; i < sSignalType.size(); ++i) {
            if (name == sSignalType[i])
                return static_cast<SignalType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
#else
    static const char* const sSignalType[4];
    static std::string signalTypeName(SignalType signalType) { return sSignalType[static_cast<int>(signalType)]; }
    static SignalType toSignalType(const std::string & name) {
        for (size_t i = 0; i < 4; ++i) {
            if (name == sSignalType[i])
                return static_cast<SignalType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
#endif

    const IndexSin::SinSource sinSource;

    const unsigned order;
    const unsigned size;
    const unsigned freq;
    const SignalType signalType;

    FFT_Signal(SignalType signalType, unsigned order, int freq, IndexSin::SinSource sinSource);
    FFT_Signal(const FFT_Signal& other);

protected:
    std::vector<VarDbl> sWave, sFreq;
};


/*
Calculated results for Signals for FFT calculation as:
    Forward:    sWave + noise -> sFrwd -> sSpec and ssSpecStep
    Roundtrip:                   sSpec -> sRound and ssRoundStep
    Reverse:    sFreq + noise -> sBack -> sRev  and ssRevStep
The Aggr signal is the average of All Sin and Cos
*/
struct FFT_Order : public FFT_Signal
{
#if __cplusplus >= 201103L
    const double NORMALIZED_ERROR_OUTLIER = 1e14;
#else
    static const double NORMALIZED_ERROR_OUTLIER;
#endif

    enum NoiseType {
        Gaussian,
        White
    };
#if __cplusplus >= 202002L
    constexpr static const std::array<std::string, 2> sNoiseType{ "Gaussian", "White" };
    constexpr static const std::string noiseTypeName(NoiseType noiseType) { return sNoiseType[static_cast<size_t>(noiseType)]; }
#if __cplusplus >= 202302L
    constexpr
#endif
    static NoiseType toNoiseType(const std::string & name) {
        for (size_t i = 0; i < sNoiseType.size(); ++i) {
            if (name == sNoiseType[i])
                return static_cast<NoiseType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
#else
    static const char* const sNoiseType[2];
    static std::string noiseTypeName(NoiseType noiseType) { return sNoiseType[static_cast<int>(noiseType)]; }
    static NoiseType toNoiseType(const std::string & name) {
        for (size_t i = 0; i < 2; ++i) {
            if (name == sNoiseType[i])
                return static_cast<NoiseType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
#endif

    enum TestType {
        Forward,
        Reverse,
        Roundtrip
    };
#if __cplusplus >= 202002L
    constexpr static const std::array<std::string, 3> sTestType{ "Forward", "Reverse", "Roundtrip" };
    constexpr static const std::string testTypeName(TestType testType) { return sTestType[static_cast<size_t>(testType)]; }
#if __cplusplus >= 202302L
    constexpr
#endif
    static TestType toTestType(const std::string & name) {
        for (size_t i = 0; i < sTestType.size(); ++i) {
            if (name == sTestType[i])
                return static_cast<TestType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
#else
    static const char* const sTestType[3];
    static std::string testTypeName(TestType testType) { return sTestType[static_cast<int>(testType)]; }
    static TestType toTestType(const std::string & name) {
        for (size_t i = 0; i < 3; ++i) {
            if (name == sTestType[i])
                return static_cast<TestType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
#endif

    struct Measure {
        Stat<double, size_t> sUncStat[3];
        Stat<double, size_t> sErrStat[3];
        Histogram<double, size_t> sHisto[3];
    };

    const NoiseType noiseType;
    const double noise;

    static std::string defaultDumpPath(unsigned minOrder=IndexSin::MIN_ORDER, unsigned maxOrder=IndexSin::MAX_ORDER);

#if __cplusplus >= 201103L
    static bool dump(std::string dumpPath = "",
                    unsigned minOrder=IndexSin::MIN_ORDER + 1, unsigned maxOrder=IndexSin::MAX_ORDER + 1,
                    std::initializer_list<IndexSin::SinSource> sSinSource = {IndexSin::SinSource::Quart, IndexSin::SinSource::Lib},
                    std::initializer_list<NoiseType> sNoiseType = {NoiseType::Gaussian, NoiseType::White},
                    std::initializer_list<double> sNoise = {0., 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.},
                    std::initializer_list<int> sFreq = {1,2,3,4,5,6,7}
                    );
#else
    static bool dump(std::string dumpPath = "",
                    unsigned minOrder=IndexSin::MIN_ORDER + 1, unsigned maxOrder=IndexSin::MAX_ORDER + 1,
                    const IndexSin::SinSource* sSinSource = 0, size_t nSinSource = 0,
                    const NoiseType* sNoiseType = 0, size_t nNoiseType = 0,
                    const double* sNoise = 0, size_t nNoise = 0,
                    const int* sFreq = 0, size_t nFreq = 0
                    );
#endif

    FFT_Order(const FFT_Signal& signal, NoiseType noiseType, double noise,
              bool traceSteps=false, unsigned minCount=64);

protected:
    std::vector<VarDbl> sFrwd, sBack;
    std::vector<VarDbl> sSpec, sRound, sRev;
    std::vector<std::vector<VarDbl> > ssSpecStep, ssRoundStep, ssRevStep;

    Measure measure;

private:
    double getNoise(Random& rand) const;
    void accum(FFT_Order::TestType testType, size_t index, const VarDbl& res, const VarDbl& err, bool hasAggr);

    static std::map<unsigned, std::map<NoiseType, std::map<double, std::map<IndexSin::SinSource, Measure> > > > ssssAggr;
    static Measure& aggr(unsigned order, NoiseType noiseType, double noise, IndexSin::SinSource sinSource);
    Measure& aggr() const { return aggr(order, noiseType, noise, sinSource); }

    void dump(std::ofstream& ofs, SignalType signalType, const Measure& measure) const;
    void dump(std::ofstream& ofs) const;
};

#if __cplusplus < 202002L
const char* const FFT_Signal::sSignalType[4] = { "Sin", "Cos", "Linear", "Aggr" };
const char* const FFT_Order::sNoiseType[2] = { "Gaussian", "White" };
const char* const FFT_Order::sTestType[3] = { "Forward", "Reverse", "Roundtrip" };
#endif
#if __cplusplus < 201103L
const double FFT_Order::NORMALIZED_ERROR_OUTLIER = 1e14;
#endif


/*
Dump the result steps for each order and sin source
*/
struct FFT_Step : public FFT_Order
{
    FFT_Step(const FFT_Signal& signal, NoiseType noiseType, double noise) :
        FFT_Order(signal, noiseType, noise, true) {}

    static std::string defaultDumpPath(IndexSin::SinSource sinSource, unsigned order);

#if __cplusplus >= 201103L
    static bool dump(IndexSin::SinSource sinSource, unsigned order,
                     std::string dumpPath = "",
                     std::initializer_list<int> sFreq = {1,2,3,4,5,6,7},
                     std::initializer_list<NoiseType> sNoiseType = {NoiseType::Gaussian},
                     std::initializer_list<double> sNoise = {0.},
                     std::string dumpOrderPath = ""
                    );
#else
    static bool dump(IndexSin::SinSource sinSource, unsigned order,
                     std::string dumpPath = "",
                     const int* sFreq = 0, size_t nFreq = 0,
                     const NoiseType* sNoiseType = 0, size_t nNoiseType = 0,
                     const double* sNoise = 0, size_t nNoise = 0,
                     std::string dumpOrderPath = ""
                    );
#endif

private:
    static void dump(std::ofstream& ofs, const std::vector<VarDbl>& sData, const std::string& context);
    void dump(std::ofstream& ofs, TestType testType,
              const std::vector<std::vector<VarDbl> >& ssStep) const;
    void dump(std::ofstream& ofs) const;
};


} // namespace var_dbl
#endif  // __Test_FFT_h__
