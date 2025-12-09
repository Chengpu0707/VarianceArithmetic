#include <map>

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
    constexpr static const std::array<std::string, 4> sSignalType{ "Sin", "Cos", "Linear", "Aggr" };
    constexpr static const std::string signalTypeName(SignalType signalType) { return sSignalType[static_cast<size_t>(signalType)]; }
    constexpr static SignalType toSignalType(const std::string & name) {
        for (size_t i = 0; i < sSignalType.size(); ++i) {
            if (name == sSignalType[i])
                return static_cast<SignalType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }
    
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
    const double NORMALIZED_ERROR_OUTLIER = 1e14;

    enum NoiseType {
        Gaussian,
        White
    };
    constexpr static const std::array<std::string, 2> sNoiseType{ "Gaussian", "White" };
    constexpr static const std::string noiseTypeName(NoiseType noiseType) { return sNoiseType[static_cast<size_t>(noiseType)]; }
    constexpr static NoiseType toNoiseType(const std::string & name) {
        for (size_t i = 0; i < sNoiseType.size(); ++i) {
            if (name == sNoiseType[i])
                return static_cast<NoiseType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }

    enum TestType {
        Forward,
        Reverse,
        Roundtrip
    };
    constexpr static const std::array<std::string, 3> sTestType{ "Forward", "Reverse", "Roundtrip" };
    constexpr static const std::string testTypeName(TestType testType) { return sTestType[static_cast<size_t>(testType)]; }
    constexpr static TestType toTestType(const std::string & name) {
        for (size_t i = 0; i < sTestType.size(); ++i) {
            if (name == sTestType[i])
                return static_cast<TestType>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SignalType";
        throw std::invalid_argument(oss.str());
    }

    struct Measure {
        Stat<double, size_t> sUncStat[3];
        Stat<double, size_t> sErrStat[3];
        Histogram<double, size_t> sHisto[3];
    };

    const NoiseType noiseType;
    const double noise;

    static std::string defaultDumpPath(unsigned minOrder=IndexSin::MIN_ORDER, unsigned maxOrder=IndexSin::MAX_ORDER);

    static bool dump(std::string dumpPath = "",
                    unsigned minOrder=IndexSin::MIN_ORDER + 1, unsigned maxOrder=IndexSin::MAX_ORDER + 1,
                    std::initializer_list<IndexSin::SinSource> sSinSource = {IndexSin::SinSource::Prec, IndexSin::SinSource::Quart, IndexSin::SinSource::Lib},
                    std::initializer_list<NoiseType> sNoiseType = {NoiseType::Gaussian, NoiseType::White}, 
                    std::initializer_list<double> sNoise = {0., 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.},
                    std::initializer_list<int> sFreq = {1,2,3,4,5,6,7} 
                    );

    FFT_Order(const FFT_Signal& signal, NoiseType noiseType, double noise, 
              bool traceSteps=false, unsigned minCount=64);

protected:
    std::vector<VarDbl> sFrwd, sBack;
    std::vector<VarDbl> sSpec, sRound, sRev;
    std::vector<std::vector<VarDbl>> ssSpecStep, ssRoundStep, ssRevStep;

    Measure measure;

private:
    double getNoise(Random& rand) const;
    void accum(FFT_Order::TestType testType, size_t index, const VarDbl& res, const VarDbl& err, bool hasAggr);

    static std::map<unsigned, std::map<NoiseType, std::map<double, std::map<IndexSin::SinSource, Measure>>>> ssssAggr;
    static Measure& aggr(unsigned order, NoiseType noiseType, double noise, IndexSin::SinSource sinSource);
    Measure& aggr() const { return aggr(order, noiseType, noise, sinSource); }

    void dump(std::ofstream& ofs, SignalType signalType, const Measure& measure) const;
    void dump(std::ofstream& ofs) const;
};


/*
Dump the result steps for each order and sin source
*/
struct FFT_Step : public FFT_Order
{
    FFT_Step(const FFT_Signal& signal, NoiseType noiseType, double noise) :
        FFT_Order(signal, noiseType, noise, true) {}

    static std::string defaultDumpPath(IndexSin::SinSource sinSource, unsigned order);

    static bool dump(IndexSin::SinSource sinSource, unsigned order,
                     std::string dumpPath = "",
                     std::initializer_list<int> sFreq = {1,2,3,4,5,6,7}, 
                     std::initializer_list<NoiseType> sNoiseType = {NoiseType::Gaussian}, 
                     std::initializer_list<double> sNoise = {0.},
                     std::string dumpOrderPath = ""
                    );

private:
    static void dump(std::ofstream& ofs, const std::vector<VarDbl>& sData, const std::string& context);
    void dump(std::ofstream& ofs, TestType testType,
              const std::vector<std::vector<VarDbl>>& ssStep) const;
    void dump(std::ofstream& ofs) const;
};







} // namespace var_dbl
#endif  // __Test_FFT_h__