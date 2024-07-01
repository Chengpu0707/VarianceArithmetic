#include "FFT.h"
#include "Test.h"
#include "ulp.h"

#include <fstream>

using namespace var_dbl;


void testBitReversion()
{
    test::assertEquals( FFT::bitReversedIndices(2), 
        std::vector{0, 2, 1, 3}, "FFT::3" );
    test::assertEquals( FFT::bitReversedIndices(3), 
        std::vector{0, 4, 2, 6, 1, 5, 3, 7}, "FFT::3");
    test::assertEquals( FFT::bitReversedIndices(4), 
        std::vector{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15}, "FFT::4" );
}

void validate(const std::vector<VarDbl>& sExpected, const std::vector<VarDbl>& sRes, double delta = 0) {
    test::assertEqual( sExpected.size(), sRes.size());
    for (int i = 0; i < sExpected.size(); ++i)
        test::assertAlmostEqual(sExpected[i].value(), sRes[i].value(), (delta > 0)? delta : sRes[i].uncertainty());
}

void testFFTOrder2Sin() 
{
    FFT fft(FFT::IndexedSin);
    const std::vector<VarDbl> sData{0,0, 1,0, 0,0, -1,0};
    const std::vector<VarDbl> sSpec{0,0, 0,2, 0,0, 0,-2}; 
    try {
        validate( sSpec, fft.transform(sData, true) );
        validate( sData, fft.transform(sSpec, false) );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}

void testFFTOrder2Cos() 
{
    FFT fft(FFT::IndexedSin);
    const std::vector<VarDbl> sData{1,0, 0,0, -1,0, 0,0};
    const std::vector<VarDbl> sSpec{0,0, 2,0, 0,0, 2,0}; 
    try {
        validate( sSpec, fft.transform(sData, true) );
        validate( sData, fft.transform(sSpec, false) );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}

void testFFTOrder3Sin() 
{ 
    const double q = std::sqrt(0.5);
    FFT fft(FFT::IndexedSin);
    const std::vector<VarDbl> sData{0.,0., q,0., 1.,0., q,0., 0,0, -q,0., -1,0, -q,0.};
    const std::vector<VarDbl> sSpec{0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4}; 
    try {
        validate( sSpec, fft.transform(sData, true), ulp(2.) );
        validate( sData, fft.transform(sSpec, false), ulp(2.) );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}

void testFFTOrder3Cos() 
{ 
    const double q = std::sqrt(0.5);
    FFT fft(FFT::IndexedSin);
    const std::vector<VarDbl> sData{1.,0., q,0., 0,0, -q,0., -1,0, -q,0., 0.,0., q,0.};
    const std::vector<VarDbl> sSpec{0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0}; 
    try {
        validate( sSpec, fft.transform(sData, true), ulp(2.) );
        validate( sData, fft.transform(sSpec, false), ulp(2.) );
    } catch (std::runtime_error e) {
        test::fail(e.what());
    }       
}


/*
 * A class to facilitate testing of FFT
 */
struct Signal {
    enum SignalType {
        Sin,
        Cos,
        Linear,   
        Aggr,  // aggregated Sin and Cos for all frequencies 
    };
    constexpr static const std::array<std::string, 4> sSignalType{ "Sin", "Cos", "Linear", "Aggr" };
    constexpr static const std::string signalTypeName(SignalType signalType) { return sSignalType[static_cast<size_t>(signalType)]; }
    
    enum NoiseType {
        Gaussian,
        White
    };
    constexpr static const std::array<std::string, 2> sNoiseType{ "Gaussian", "White" };
    constexpr static const std::string noiseTypeName(NoiseType noiseType) { return sNoiseType[static_cast<size_t>(noiseType)]; }

    enum TestType {
        Forward,
        Reverse,
        Roundtrip
    };
    constexpr static const std::array<std::string, 3> sTestType{ "Forward", "Reverse", "Roundtrip" };
    constexpr static const std::string testTypeName(TestType testType) { return sTestType[static_cast<size_t>(testType)]; }

    struct Measure {
        Stat<double> sStat[3];
        Histogram sHisto[3];
    };

    static std::map<unsigned, std::map<NoiseType, std::map<double, std::map<FFT::SinSource, Measure>>>> ssssAggr;
    static void dump();

    const FFT::SinSource sinSource;
    const FFT fft;

    const unsigned order;
    const unsigned size;
    const unsigned freq;
    const SignalType signalType;
    std::vector<VarDbl> sWave, sFreq;

    const NoiseType noiseType;
    const double noise;

    std::vector<VarDbl> sData, sSpec, sRound, sBack, sRev;
    Measure measure;

    Signal(SignalType signalType, unsigned order, int freq, NoiseType noiseType, double noise, FFT::SinSource sinSource);

private:
    static Measure& aggr(unsigned order, NoiseType noiseType, double noise, FFT::SinSource sinSource);
    Measure& aggr() const { return aggr(order, noiseType, noise, sinSource); }
    double getNoise(Random& rand) const;

    static void run(std::ofstream& of, FFT::SinSource sinSource, NoiseType noiseType, double noise, 
            SignalType signalType, int order, int freq);
    static void dump(std::ofstream& of, FFT::SinSource sinSource, NoiseType noiseType, double noise, 
            SignalType signalType, int order, int freq, const Measure& measure);

};


std::map<unsigned, std::map<Signal::NoiseType, std::map<double, std::map<FFT::SinSource, Signal::Measure>>>> Signal::ssssAggr;


double Signal::getNoise(Random& rand) const
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


Signal::Measure& Signal::aggr(unsigned order, NoiseType noiseType, double noise, FFT::SinSource sinSource)
{
    auto it1 = ssssAggr.find(order);
    if (it1 == ssssAggr.end()) {
        auto res = ssssAggr.insert({order, 
            std::map<NoiseType, std::map<double, std::map<FFT::SinSource, Measure>>>()});
        it1 = res.first;
    }
    auto it2 = it1->second.find(noiseType);
    if (it2 == it1->second.end()) {
        auto res = it1->second.insert({noiseType, std::map<double, std::map<FFT::SinSource, Measure>>()});
        it2 = res.first;
    }
    auto it3 = it2->second.find(noise);
    if (it3 == it2->second.end()) {
        auto res = it2->second.insert({noise, std::map<FFT::SinSource, Measure>()});
        it3 = res.first;
    }
    auto it4 = it3->second.find(sinSource);
    if (it4 == it3->second.end()) {
        auto res = it3->second.insert({sinSource, Measure()});
        it4 = res.first;
    }
    return it4->second;
}


Signal::Signal(SignalType signalType, unsigned order, int freq, NoiseType noiseType, double noise, FFT::SinSource sinSource) :
	order(order), freq(freq), signalType(signalType), noiseType(noiseType), noise(std::abs(noise)), sinSource(sinSource),
	fft(sinSource), size(1 << order)
{
	std::ostringstream oss;
	if (freq * 2 >= size) {
	    oss << "Sin invalid freq " << freq << " for order " << order;
	    throw std::invalid_argument(oss.str());
	}
	switch (noiseType) {
	case Gaussian:
		break;
	case White:
		break;
	default:
		oss << "Unknown noise type " << noiseType;
		throw std::invalid_argument(oss.str());
	}
	Random rand(0, noise);      

	sWave.reserve(size << 1);
	const int peak = size >> 1;
	switch (signalType) {
	case Sin:
		for (long i = 0; i < size; ++i) {
		    sWave.push_back(fft.sin(freq * i, order));
		    sWave.push_back( 0 );
		}
		sFreq.insert(sFreq.end(), size << 1, 0);
		sFreq[(freq << 1) + 1] = peak;
		sFreq[((size - freq) << 1) + 1] = - peak;
		break;
	case Cos:
		for (long i = 0; i < size; ++i) {
		    sWave.push_back(fft.cos(freq * i, order));
		    sWave.push_back( 0 );
		}
		sFreq.insert(sFreq.end(), size << 1, 0);
		sFreq[freq << 1] = peak;
		sFreq[(size - freq) << 1] = peak;
		break;
	case Linear:
		if (order == FFT::MAX_ORDER) {
		    oss << "FFT order " << order << " is too large for FFT sine order " << FFT::MAX_ORDER;
		    throw std::invalid_argument(oss.str());
		}
		sFreq.reserve(size << 1);
		for (int i = 0; i < size; ++i) {
		    sWave[i << 1] = i;
		    sWave[(i << 1) + 1] = 0;
		    if (i == 0) {
                sFreq.push_back(size * (size - 1) / 2);
                sFreq.push_back(0);
		    } else {
                sFreq[i << 1] = -peak;
                sFreq[(i << 1) + 1] = -peak / fft.sin(i, order + 1) * fft.cos(i, order + 1);
		    }
	    }
		break;
	default:
		oss << "FFT signal Type " << signalType << " is unknown ";
		throw std::invalid_argument(oss.str());
	}

	sData.reserve(size << 1);
	sBack.reserve(size << 1);
	for (int i = 0; i < (size << 1); ++i) { 
        sData.push_back(VarDbl(sWave[i].value() + getNoise(rand), noise*noise + sWave[i].variance(), true));
        sBack.push_back(VarDbl(sFreq[i].value() + getNoise(rand), noise*noise + sFreq[i].variance(), true));
	}

	sSpec = fft.transform(sData, true);
	sRound = fft.transform(sSpec, false);
	sRev = fft.transform(sBack, false);
    bool hasAggr = (signalType == Sin) || (signalType == Cos);
                                                                                                                                                               
	for (int i = 0; i < (size << 1); ++i) {   
	    const double unc1 = sSpec[i].uncertainty();
	    measure.sStat[Forward].add(unc1);
	    if (hasAggr)
		    aggr().sStat[Forward].add(unc1);
	    if (unc1 > 0) {
		    measure.sHisto[Forward].add((sSpec[i].value() - sFreq[i].value())/unc1);
            if (hasAggr)
                aggr().sHisto[Forward].add((sSpec[i].value() - sFreq[i].value())/unc1);
	    }

	    const double unc2 = sRound[i].uncertainty();
	    measure.sStat[Roundtrip].add(unc2);
	    if (hasAggr)
		    aggr().sStat[Roundtrip].add(unc2);
	    if (unc2 > 0) {
		    measure.sHisto[Roundtrip].add((sRound[i].value() - sData[i].value())/unc2);
            if (hasAggr)
                aggr().sHisto[Roundtrip].add((sRound[i].value() - sData[i].value())/unc2);
	    }

	    const double unc3 = sRev[i].uncertainty();
	    measure.sStat[Reverse].add(unc3);
	    if (hasAggr)
		    aggr().sStat[Reverse].add(unc3);
	    if (unc3 > 0) {
            measure.sHisto[Reverse].add((sRev[i].value() - sWave[i].value())/unc3);
            if (hasAggr)
                aggr().sHisto[Reverse].add((sRev[i].value() - sWave[i].value())/unc3);
	    }
	}
}


void Signal::dump(std::ofstream& of, FFT::SinSource sinSource, NoiseType noiseType, double noise, SignalType signalType, int order, int freq,
        const Measure& measure)
{
    for (unsigned testType = 0; testType < Signal::sTestType.size(); ++testType) {
        const Stat<double>& stat( measure.sStat[testType] );
        const Histogram& histo( measure.sHisto[testType] );
        of << FFT::sinSourceName(sinSource) << "\t" << Signal::noiseTypeName(noiseType) << "\t" << noise
            << "\t" << Signal::signalTypeName(signalType) << "\t" << order << "\t" << freq << "\t" << Signal::sTestType[testType]
            << "\t" << stat.mean() << "\t" << stat.std() << "\t" << stat.min() << "\t" << stat.max() 
            << "\t" << histo.mean() << "\t" << histo.std() << "\t" << histo.min() << "\t" << histo.max() 
            << histo.formatted() << "\n";
    }
}


void Signal::run(std::ofstream& of, FFT::SinSource sinSource, NoiseType noiseType, double noise, SignalType signalType, int order, int freq)
{
    const Signal signal(signalType, order, freq, noiseType, noise, sinSource);
    dump(of, sinSource, noiseType, noise, signalType, order, freq, signal.measure);
}


void Signal::dump() 
{
    Histogram histo;
    std::ofstream of("./Output/FFT_4_18.txt");
    of << "SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest"
        << "\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Maximum"
        << "\tError Mean\tError Deviation\tError Minimum\tError Maximum"
        << histo.header() << "\n";
        
    for (int order = 4; order < 19; ++order) {
        for (int i = -17; i <= 0; ++i) {
            const double noise = std::pow(10, i);
            for (unsigned noiseType = 0; noiseType < sNoiseType.size(); ++noiseType) {
                for (unsigned sinSource = 0; sinSource < FFT::sSinSource.size(); ++sinSource) {
                    run(of, FFT::SinSource(sinSource), NoiseType(noiseType), noise, Signal::Linear, order, 0);
                    for (int freq = 1; freq < 8; ++freq) {
                        run(of, FFT::SinSource(sinSource), NoiseType(noiseType), noise, Signal::Sin, order, freq);
                        run(of, FFT::SinSource(sinSource), NoiseType(noiseType), noise, Signal::Cos, order, freq);
                    }
                    dump(of, FFT::SinSource(sinSource), NoiseType(noiseType), noise, Signal::Aggr, order, 0, 
                         aggr(order, NoiseType(noiseType), noise, FFT::SinSource(sinSource)));  
                }
                of.flush();
            }
        }
    }
}


int main()
{
    testBitReversion();

    testFFTOrder2Sin();
    testFFTOrder2Cos();
    testFFTOrder3Sin();
    testFFTOrder3Cos();

    Signal::dump();

    std::cout << "All FFT tests are successful";
    return 0;
}


