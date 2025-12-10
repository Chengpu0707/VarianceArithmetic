#include "Test.h"
#include "TestFFT.h"
#include "ulp.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace var_dbl;

void validate(const std::vector<VarDbl>& sExpected, const std::vector<VarDbl>& sRes, double delta = 0) {
    test::assertEqual( sExpected.size(), sRes.size());
    for (int i = 0; i < sExpected.size(); ++i)
        test::assertAlmostEqual(sExpected[i].value(), sRes[i].value(), (delta > 0)? delta : sRes[i].uncertainty());
}

void testFFTOrder2Sin(const IndexSin::SinSource src, double delta) 
{
    FFT fft(src);
    const std::vector<VarDbl> sFrwd{0,0, 1,0, 0,0, -1,0};
    const std::vector<VarDbl> sSpec{0,0, 0,2, 0,0, 0,-2}; 
    try {
        validate( sSpec, fft.transform(sFrwd, true, true), delta);
        validate( sFrwd, fft.transform(sSpec, false, true), delta );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}

void testFFTOrder2Cos(const IndexSin::SinSource src, double delta) 
{
    FFT fft(src);
    const std::vector<VarDbl> sFrwd{1,0, 0,0, -1,0, 0,0};
    const std::vector<VarDbl> sSpec{0,0, 2,0, 0,0, 2,0}; 
    try {
        validate( sSpec, fft.transform(sFrwd, true, true),  delta );
        validate( sFrwd, fft.transform(sSpec, false, true), delta );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}

void testFFTOrder3Sin(const IndexSin::SinSource src, double delta) 
{ 
    const double q = std::sqrt(0.5);
    FFT fft(src);
    const std::vector<VarDbl> sFrwd{0.,0., q,0., 1.,0., q,0., 0,0, -q,0., -1,0, -q,0.};
    const std::vector<VarDbl> sSpec{0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4}; 
    try {
        validate( sSpec, fft.transform(sFrwd, true, true),  delta );
        validate( sFrwd, fft.transform(sSpec, false, true), delta );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}

void testFFTOrder3Cos(const IndexSin::SinSource src, double delta) 
{ 
    const double q = std::sqrt(0.5);
    FFT fft(src);
    const std::vector<VarDbl> sFrwd{1.,0., q,0., 0,0, -q,0., -1,0, -q,0., 0.,0., q,0.};
    const std::vector<VarDbl> sSpec{0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0}; 
    try {
        validate( sSpec, fft.transform(sFrwd, true, true),  delta );
        validate( sFrwd, fft.transform(sSpec, false, true), delta );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}


void testBitReversion()
{
    test::assertEquals( FFT::bitReversedIndices(2), 
        std::vector{0, 2, 1, 3}, "FFT::3" );
    test::assertEquals( FFT::bitReversedIndices(3), 
        std::vector{0, 4, 2, 6, 1, 5, 3, 7}, "FFT::3");
    test::assertEquals( FFT::bitReversedIndices(4), 
        std::vector{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15}, "FFT::4" );
}



/*
 * A class to facilitate testing of FFT
 */

FFT_Signal::FFT_Signal(SignalType signalType, unsigned order, int freq, IndexSin::SinSource sinSource) :
	order(order), size(1 << order), freq(freq), signalType(signalType), sinSource(sinSource), FFT(sinSource)
{
	std::ostringstream oss;
	if (freq * 2 >= size) {
	    oss << "Sin invalid freq " << freq << " for order " << order;
	    throw std::invalid_argument(oss.str());
	}

	sWave.reserve(size << 1);
	const long peak = size >> 1;
	switch (signalType) {
	case Sin:
		for (long i = 0; i < (size << 1); i += 2) {
		    sWave.push_back(_sin.sin(freq * i, order));
		    sWave.push_back(0);
		}
		sFreq.insert(sFreq.end(), size << 1, 0);
		sFreq[(freq << 1) + 1] = peak;
		sFreq[((size - freq) << 1) + 1] = - peak;
		break;
	case Cos:
		for (long i = 0; i < (size << 1); i += 2) {
		    sWave.push_back(_sin.cos(freq * i, order));
		    sWave.push_back(0);
		}
		sFreq.insert(sFreq.end(), size << 1, 0);
		sFreq[freq << 1] = peak;
		sFreq[(size - freq) << 1] = peak;
		break;
	case Linear:
		sFreq.reserve(size << 1);
		for (long i = 0; i < size; ++i) {
		    sWave.push_back(i);
            sWave.push_back(0);
		    if (i == 0) {
                sFreq.push_back(((double) size) * (size - 1) / 2);
                sFreq.push_back(0);
		    } else {
                sFreq.push_back(-peak);
                sFreq.push_back(-peak * _sin.cos(i, order) / _sin.sin(i, order));
		    }
	    }
		break;
	default:
		oss << "FFT signalType Type " << signalType << " is unknown ";
		throw std::invalid_argument(oss.str());
	}
}


FFT_Signal::FFT_Signal(const FFT_Signal& other) :
    FFT(other.sinSource),
    sinSource(other.sinSource), order(other.order), size(other.size), freq(other.freq), signalType(other.signalType),
    sWave(other.sWave), sFreq(other.sFreq)
{
}


std::map<unsigned, std::map<FFT_Order::NoiseType, std::map<double, std::map<IndexSin::SinSource, FFT_Order::Measure>>>> FFT_Order::ssssAggr;


FFT_Order::Measure& FFT_Order::aggr(unsigned order, NoiseType noiseType, double noise, IndexSin::SinSource sinSource)
{
    auto it1 = ssssAggr.find(order);
    if (it1 == ssssAggr.end()) {
        auto res = ssssAggr.insert({order, 
            std::map<NoiseType, std::map<double, std::map<IndexSin::SinSource, Measure>>>()});
        it1 = res.first;
    }
    auto it2 = it1->second.find(noiseType);
    if (it2 == it1->second.end()) {
        auto res = it1->second.insert({noiseType, std::map<double, std::map<IndexSin::SinSource, Measure>>()});
        it2 = res.first;
    }
    auto it3 = it2->second.find(noise);
    if (it3 == it2->second.end()) {
        auto res = it2->second.insert({noise, std::map<IndexSin::SinSource, Measure>()});
        it3 = res.first;
    }
    auto it4 = it3->second.find(sinSource);
    if (it4 == it3->second.end()) {
        auto res = it3->second.insert({sinSource, Measure()});
        it4 = res.first;
    }
    return it4->second;
}


double FFT_Order::getNoise(Random& rand) const
{
    switch (noiseType) {
        case Gaussian:
            return rand.gauss();
        case White:
            return rand.white();
        default:
            return 0;
    }
}

void FFT_Order::accum(TestType testType, size_t index, const VarDbl& res, const VarDbl& err, bool hasAggr)
{
    measure.sUncStat[testType].addAt(res.uncertainty(), index);
    measure.sErrStat[testType].addAt(err.value(), index);
    if (hasAggr) {
        aggr().sUncStat[testType].addAt(res.uncertainty(), index);
        aggr().sErrStat[testType].addAt(err.value(), index);
    }
    if (err.uncertainty() > 0) {
        const double norm = err.value() / err.uncertainty();
        if (abs(norm) < FFT_Order::NORMALIZED_ERROR_OUTLIER) {
            measure.sHisto[testType].addAt(norm, index);
            if (hasAggr) {
                aggr().sHisto[testType].addAt(norm, index);
            }
        } else
            std::cerr << "For signal=" << signalTypeName(signalType) << " order=" << order << " freq=" << freq 
                      << " NoiseType=" << noiseTypeName(noiseType) << " noise=" << noise << " test=" << testTypeName(testType)
                      << " index=" << index << ", normaliized error outlier " << norm << " for " << err << " from " << res << std::endl;
    }
}
  

FFT_Order::FFT_Order(const FFT_Signal& signal, NoiseType noiseType, double noise, 
                     bool traceSteps, unsigned minCount) :
        FFT_Signal(signal), 
        noiseType(noiseType), noise(std::abs(noise))  
{
	switch (noiseType) {
	case Gaussian:
		break;
	case White:
		break;
	default:
        std::ostringstream oss;
		oss << "Unknown noise type " << noiseType;
		throw std::invalid_argument(oss.str());
	}
	Random rand(0, noise);      

	sFrwd.reserve(size << 1);
	sBack.reserve(size << 1);
    do {
        sFrwd.clear();
        sBack.clear();
        for (int i = 0; i < (size << 1); ++i) { 
            sFrwd.push_back(sWave[i] + VarDbl(getNoise(rand), noise));
            sBack.push_back(sFreq[i] + VarDbl(getNoise(rand), noise));
        }

        sSpec = transform(sFrwd, true, traceSteps);
        ssSpecStep = ssStep;
        sRound = transform(sSpec, false, traceSteps);
        ssRoundStep = ssStep;
        sRev = transform(sBack, false, traceSteps);
        ssRevStep = ssStep;

        if (traceSteps) {
            ssSpecStep.push_back(sFreq);
            ssSpecStep.push_back(std::vector<VarDbl>());
            ssRoundStep.push_back(sFrwd);
            ssRoundStep.push_back(std::vector<VarDbl>());
            ssRevStep.push_back(sWave);
            ssRevStep.push_back(std::vector<VarDbl>());
        }
        const bool hasAggr = (signalType == Sin) || (signalType == Cos);
        for (size_t i = 0; i < (size << 1); ++i) { 
            const VarDbl errSpec = sSpec[i] - sFreq[i];
            accum(Forward, i, sSpec[i], errSpec, hasAggr);
            if (traceSteps)
                ssSpecStep.back().push_back(errSpec);
            const VarDbl errRound = sRound[i] - sFrwd[i];
            accum(Roundtrip, i, sRound[i], errRound, hasAggr);
            if (traceSteps)
                ssRoundStep.back().push_back(errRound);
            const VarDbl errRev = sRev[i] - sWave[i];
            accum(Reverse, i, sRev[i], errRev, hasAggr);
            if (traceSteps)
                ssRevStep.back().push_back(errRev);
        }
    } while (!traceSteps && (0 < noise) && (measure.sUncStat[Roundtrip].count() < minCount));
}


void FFT_Order::dump(std::ofstream& ofs, SignalType signalType, const Measure& measure) const
{
    for (unsigned testType = 0; testType < FFT_Order::sTestType.size(); ++testType) {
        const Stat<double, size_t>& uncStat( measure.sUncStat[testType] );
        const Stat<double, size_t>& errStat( measure.sUncStat[testType] );
        const Histogram<double, size_t>& histo( measure.sHisto[testType] );
        ofs << IndexSin::sinSourceName(sinSource) << '\t' << FFT_Order::noiseTypeName(noiseType) << '\t' << noise
            << '\t' << FFT_Signal::signalTypeName(signalType) << '\t' << order << '\t' << freq << '\t' << FFT_Order::sTestType[testType]
            << '\t' << uncStat.count()<< '\t' << uncStat.mean() << '\t' << uncStat.std() 
                << '\t' << uncStat.min() << '\t' << (uncStat.minAt().has_value()? std::to_string(uncStat.minAt().value()) : "") 
                << '\t' << uncStat.max() << '\t' << (uncStat.maxAt().has_value()? std::to_string(uncStat.maxAt().value()) : "")  
            << '\t' << errStat.count()<< '\t' << errStat.mean() << '\t' << errStat.std() 
                << '\t' << errStat.min() << '\t' << (errStat.minAt().has_value()? std::to_string(errStat.minAt().value()) : "") 
                << '\t' << errStat.max() << '\t' << (errStat.maxAt().has_value()? std::to_string(errStat.maxAt().value()) : "")  
            << '\t' << histo.count() << '\t' << histo.mean() << '\t' << histo.std() 
                << '\t' << histo.min() << '\t' << (histo.minAt().has_value()? std::to_string(histo.minAt().value()) : "")  
                << '\t' << histo.max() << '\t' << (histo.maxAt().has_value()? std::to_string(histo.maxAt().value()) : "")
            << '\t' << histo.lowers() << '\t' << histo.uppers()
            << histo.formatted() << "\n";
    }
}

void FFT_Order::dump(std::ofstream& ofs) const
{
    dump(ofs, signalType, measure);
}

std::string FFT_Order::defaultDumpPath(unsigned minOrder, unsigned maxOrder)
{
    std::ostringstream oss;
    oss << "./Output/FFT_" << minOrder << "_"<< maxOrder <<".txt";
    return oss.str();
}

bool FFT_Order::dump(
            std::string dumpPath,
            unsigned minOrder, unsigned maxOrder, 
            std::initializer_list<IndexSin::SinSource> sSinSource,
            std::initializer_list<NoiseType> sNoiseType, 
            std::initializer_list<double> sNoise,
            std::initializer_list<int> sFreq 
    ) 
{
    Histogram histo;
    std::ostringstream oss;
    oss << "SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest"
        << "\tUncertainty Count\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Minimum At\tUncertainty Maximum\tUncertainty Maximum At"
        << "\tValue Count\tValue Mean\tValue Deviation\tValue Minimum\tValue Minimum At\tValue Maximum\tValue Maximum At"
        << "\tError Count\tError Mean\tError Deviation\tError Minimum\tError Minimum At\tError Maximum\tError Maximum At"
        << "\tLower Count\tUpper Count" << histo.header();
    const std::string HEADER = oss.str();
    if (dumpPath.empty())
        dumpPath = defaultDumpPath(minOrder, maxOrder);

    std::unordered_map<int, std::unordered_map<NoiseType, std::unordered_map<double, 
            std::unordered_map<IndexSin::SinSource, std::unordered_set<TestType>>>>> sssssAggr;
    std::ifstream ifs(dumpPath);
    std::string line;
    size_t n;
    if (ifs.is_open() && std::getline(ifs, line) && (line == HEADER)) {
        for (n = 1; std::getline(ifs, line); ++n) {
            std::istringstream iss(line);
            std::string sinName, noiseName, signalName, testName;
            int order, freq;
            double noise;
            iss >> sinName >> noiseName >> noise >> signalName >> order >> freq >> testName;
            if ((noise < 0) || (1 < noise)) {
                std::cout << "At line #" << n <<", invalid noise=" << noise << " in: " << line << std::endl;
                return false;
            }
            if ((order < IndexSin::MIN_ORDER) || (IndexSin::MAX_ORDER < order)) {
                std::cout << "At line #" << n <<", invalid order=" << order << " in: " << line << std::endl;
                return false;
            }
            if ((freq < 0) || (7 < freq)) {
                std::cout << "At line #" << n <<", invalid freq=" << freq << " in: " << line << std::endl;
                return false;
            }
            try {
                IndexSin::SinSource sinSource = IndexSin::toSinSource(sinName);
                NoiseType noiseType = FFT_Order::toNoiseType(noiseName);
                SignalType signalType = FFT_Signal::toSignalType(signalName);
                TestType testType = FFT_Order::toTestType(testName);
                if (std::string(FFT_Signal::sSignalType[signalType]) != "Aggr")
                    continue;
                auto itO = sssssAggr.find(order);
                if (itO == sssssAggr.end()) {
                    auto res = sssssAggr.insert({order, std::unordered_map<NoiseType, std::unordered_map<double, 
                                        std::unordered_map<IndexSin::SinSource, std::unordered_set<TestType>>>>()});
                    itO = res.first;
                }
                auto itN = itO->second.find(noiseType);
                if (itN == itO->second.end()) {
                    auto res = itO->second.insert({noiseType, std::unordered_map<double, 
                                        std::unordered_map<IndexSin::SinSource, std::unordered_set<TestType>>>()});
                    itN = res.first;
                }
                auto itn = itN->second.find(noise);
                if (itn == itN->second.end()) {
                    auto res = itN->second.insert({noise, std::unordered_map<IndexSin::SinSource, std::unordered_set<TestType>>()});
                    itn = res.first;
                }
                auto itS = itn->second.find(sinSource);
                if (itS == itn->second.end()) {
                    auto res = itn->second.insert({sinSource, std::unordered_set<TestType>()});
                    itS = res.first;
                }
                auto res = itS->second.insert(testType);
                if (!res.second) {
                    std::cout << "At line #" << n << ", duplicated order=" << order << ", noiseType=" << noiseType << ", noise=" << noise 
                            << ", sinSource=" << sinSource << ", signalType=" << signalType
                            << ": " << line << std::endl;
                    return false;
                }
            } catch (std::invalid_argument ex) {
                std::cout  << "At line #" << n << " " << ex.what() << ": " << line << std::endl;
                return false;
            }
        }
    }

    std::ofstream ofs(dumpPath, sssssAggr.empty()? std::ios::out : std::ios::app);
    if (!ofs) {
        std::cout << "Fail to open " << dumpPath;
        return false;
    }
    ofs << std::scientific << std::setprecision(20);
    if (sssssAggr.empty())
        ofs << HEADER << "\n";   
        
    for (NoiseType noiseType: sNoiseType) {
        for (unsigned order = minOrder; order < maxOrder; ++order) {
            const unsigned half = 1 << (order - 1);
            for (IndexSin::SinSource sinSource: sSinSource) {
                std::vector<FFT_Signal> sSignal;
                for (double noise: sNoise) {
                    if (auto itO = sssssAggr.find(order); itO != sssssAggr.end()) {
                        if (auto itN = itO->second.find(noiseType); itN != itO->second.end()) {
                            if (auto itn = itN->second.find(noise); itn != itN->second.end()) {                               
                                if (auto itS = itn->second.find(sinSource); itS != itn->second.end()) {
                                    if (itS->second.size() == 3)
                                        continue;
                                }
                            }
                        }
                    }
                    const std::time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                    std::cout << std::put_time(std::localtime(&currentTime), "%Y-%m-%d %H:%M:%S") << ": "
                            << "Start calulation order=" << order << ", sinSource=" << IndexSin::sSinSource[sinSource] 
                            << ", noiseType=" << FFT_Order::sNoiseType[noiseType] <<", noise=" << noise << std::endl;
                    if (sSignal.empty()) {
                        for (int freq : sFreq) {
                            if (freq >= half)
                                break;
                            sSignal.emplace_back(FFT_Signal::Sin, order, freq, sinSource);
                            sSignal.emplace_back(FFT_Signal::Cos, order, freq, sinSource);
                        }
                        sSignal.emplace_back(FFT_Signal::Linear, order, 0, sinSource);
                        const std::time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                        std::cout << std::put_time(std::localtime(&currentTime), "%Y-%m-%d %H:%M:%S") << ": "
                                  << "Finish create signal for order=" << order << ", sinSource=" << IndexSin::sSinSource[sinSource] << std::endl;
                    }
                    for (const FFT_Signal& signal: sSignal) {
                        FFT_Order calc(signal, NoiseType(noiseType), noise);
                        calc.dump(ofs);
                        if (signal.signalType == SignalType::Linear)
                            calc.dump(ofs, SignalType::Aggr, calc.aggr());
                        ofs.flush();
                    }                   
                }
            }
        }
    }
    return true;
}

void FFT_Step::dump(std::ofstream& ofs, const std::vector<VarDbl>& sData, const std::string& context)
{
    for (int imag: {0, 1}) {
        for (int value: {1, 0}) {
            ofs << context << '\t' << imag << '\t' << value;;
            for (size_t i = imag; i < sData.size(); i += 2) {
                if (value)
                    ofs << '\t' << sData[i].value();
                else
                    ofs << '\t' << sData[i].uncertainty();
            }
            ofs << '\n';
        }
    }
    ofs.flush();
}


void FFT_Step::dump(std::ofstream& ofs, TestType testType,
        const std::vector<std::vector<VarDbl>>& ssStep) const
{
    std::ostringstream oss;
    for (size_t step = 0; step < ssStep.size(); ++step) {
        oss.str("");
        oss << IndexSin::sinSourceName(sinSource) << '\t' << noiseTypeName(noiseType) << '\t' << noise
            << '\t' << FFT_Signal::signalTypeName(signalType) << '\t' << order << '\t' << freq << '\t' 
            << FFT_Order::testTypeName(testType) << '\t' << step;
        dump(ofs, ssStep[step], oss.str());
    }
    const std::vector<VarDbl>& sCompare((testType == TestType::Forward)? sFreq : (testType == TestType::Reverse)? sWave : sFrwd);
}


void FFT_Step::dump(std::ofstream& ofs) const
{
    dump(ofs, TestType::Forward, ssSpecStep);
    dump(ofs, TestType::Roundtrip, ssRoundStep);
    dump(ofs, TestType::Reverse, ssRevStep);
}


std::string FFT_Step::defaultDumpPath(IndexSin::SinSource sinSource, unsigned order)
{
    std::ostringstream oss;
    oss << "./Output/FFT_Step_" << order << "_" << IndexSin::sinSourceName(sinSource) << ".txt";
    return oss.str();
}

bool FFT_Step::dump(
        IndexSin::SinSource sinSource, unsigned order, std::string dumpPath, 
        std::initializer_list<int> sFreq, 
        std::initializer_list<NoiseType> sNoiseType, std::initializer_list<double> sNoise,
        std::string dumpOrderPath
    )
{
    if (order < IndexSin::MIN_ORDER || order > IndexSin::MAX_ORDER) {
        std::cerr << "Invalid order=" << order << " for SinSource=" << IndexSin::sinSourceName(sinSource);
        return false;
    }
    std::cout << "Start calculate FFT step order=" << order << " freq=[";
    for (int freq: sFreq) 
        std::cout << freq << ',';
    std::cout << "] sinSource=" << IndexSin::sinSourceName(sinSource) << " noise=[";
    for (double noise : sNoise)
        std::cout << noise << ','; 
    std::cout << "] with order result to " << dumpOrderPath << std::endl;
    const size_t size = 1 << order;

    std::ostringstream oss;
    if (dumpPath.empty())
        dumpPath = defaultDumpPath(sinSource, order);
    std::ofstream ofs(dumpPath);
    if (!ofs.is_open()) {
        std::cerr << "Failed to open file " << dumpPath << " to dump FFT with SinSource=" << IndexSin::sinSourceName(sinSource);
        return false;
    }
    ofs << std::scientific << std::setprecision(20);
    ofs << "SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest\tStep\tImag\tValue";
    for (size_t i = 0; i < size; ++i) {
        ofs << '\t' << i;
    }
    ofs << '\n';
    std::vector<VarDbl> sCosSin;
    sCosSin.reserve(size);
    const IndexSin sin(sinSource);
    for (size_t i = 0; i < size; ++i) {
        sCosSin.push_back(sin.cos(i, order));
        sCosSin.push_back(sin.sin(i, order));
    }
    oss.str("");
    oss << IndexSin::sinSourceName(sinSource) << "\t\t0\t\t" << order << "\t\tCosSin\t";
    dump(ofs, sCosSin, oss.str());

    std::vector<FFT_Signal> sSignal;
    sSignal.reserve(sFreq.size() * 2 + 1);
    for (const SignalType signalType : {SignalType::Sin, SignalType::Cos}) {
        for (int freq : sFreq) {
            if (freq >= (1 << (order - 1)))
                continue;
            sSignal.emplace_back(signalType, order, freq, sinSource);
        }
    }
    sSignal.emplace_back(SignalType::Linear, order, 0, sinSource);

    for (NoiseType noiseType : sNoiseType) {
        for (double noise : sNoise) {
            for (const FFT_Signal& signal : sSignal) {
                FFT_Step calc(signal, noiseType, noise);
                calc.dump(ofs);
            }
        }
    }
    if (dumpOrderPath.empty())
        return true;
    std::cout << "Start dump FFT order to " << dumpOrderPath << std::endl;
    return FFT_Order::dump(dumpOrderPath, order, order + 1, {sinSource}, sNoiseType, sNoise, sFreq);
}

struct FFT_Test: public FFT_Order {
    FFT_Test(const FFT_Signal& signal, NoiseType noiseType, double noise) :
        FFT_Order(signal, noiseType, noise, true)
    {}

    IndexSin getIndexSin() const { return _sin; }
    std::vector<std::vector<VarDbl>> getSpecStep() const {return ssSpecStep; }
};


int main(int argc, char* argv[])
{
    for (auto name: IndexSin::sSinSource) {
        IndexSin::SinSource src = IndexSin::toSinSource(name);
        testFFTOrder2Sin(src, ulp(2.0));
        testFFTOrder2Cos(src, ulp(2.0));
        testFFTOrder3Sin(src, ulp(4.0));
        testFFTOrder3Cos(src, ulp(4.0));
    }

    testBitReversion();

    if (argc == 1) {
        std::cout << "Start calculate FFT step" << std::endl;
        for (size_t order = 2; order < 7; ++order) {
            std::cout << "Start calculate FFT step order=" << order << std::endl;
            for (auto src : {"Prec", "Quart", "Lib"}) {
                FFT_Step::dump(IndexSin::toSinSource(src), order);
            }
        }

        std::cout << "Start calculate FFT order" << std::endl;
        const std::string dumpPath = FFT_Order::defaultDumpPath(2, 6);
        std::remove(dumpPath.c_str());
        test::assertTrue(FFT_Order::dump(dumpPath, 2, 6));
        std::cout << "Start reading FFT order" << std::endl;
        test::assertTrue(FFT_Order::dump("", 2, 6));
    } else if ((argc == 2) && (std::string(argv[1]) == "Test"))
        test::assertTrue(FFT_Order::dump());
    else if (argc <= 5) {
        const unsigned order = atoi(argv[1]);
        const double noise = (2 < argc)? atof(argv[2]) : 1e-3;
        IndexSin::SinSource sinSource = (3 < argc)? IndexSin::toSinSource(argv[3]) : IndexSin::SinSource::Quart;
        const int freq = (4 < argc)? atoi(argv[4]): 1;
        std::ostringstream oss;
        oss << "./Output/FFT_Step_" << order << '_' << freq << '_' 
            << IndexSin::sinSourceName(sinSource) << '_' << noise << ".txt";
        const std::string dumpStepPath = oss.str();
        oss.str("");
        oss << "./Output/FFT_Order_" << order << '_' << freq << '_' 
            << IndexSin::sinSourceName(sinSource) << '_' << noise << ".txt";
        const std::string dumpOrderPath = oss.str();
        FFT_Step::dump(sinSource, order,
                    dumpStepPath, {freq}, {FFT_Order::NoiseType::Gaussian}, {noise},
                    dumpOrderPath);
    } else {
        std::cout << "Invalid argument:";
        for (unsigned i = 1; i < argc; ++i) {
            std::cout << '\t' << argv[i];
        }
        std::cout << '\n';
    }

    std::cout << "All FFT tests are successful" << std::endl;
    return 0;
}


