/*
A sin() with index frequence as input, with resolution \pi/size()
size() = 1 << "order".
When the index "freq" == size, the result is sin(\pi).
This function has minimal float error when "withUncertainty"==False.

When "withUncertainty"==True, use regression to calculate the sin

*/
#include "Stat.h"
#include "VarDbl.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#if __cplusplus >= 202002L
#include <numbers>
#endif
#if __cplusplus >= 201103L
#include <array>
#include <string>
#endif
#include <sstream>
#include <vector>


#ifndef __IndexSin_h__
#define __IndexSin_h__
namespace var_dbl
{

#if __cplusplus < 202002L
namespace _detail {
#if __cplusplus >= 201103L
    constexpr double _pi    = 3.14159265358979323846;
    constexpr long double _pi_ld = 3.14159265358979323846264338327950288L;
#else
    static const double _pi    = 3.14159265358979323846;
    static const long double _pi_ld = 3.14159265358979323846264338327950288L;
#endif
} // namespace _detail
#define _VAR_DBL_PI    var_dbl::_detail::_pi
#define _VAR_DBL_PI_LD var_dbl::_detail::_pi_ld
#else
#define _VAR_DBL_PI    std::numbers::pi
#define _VAR_DBL_PI_LD std::numbers::pi_v<long double>
#endif


class IndexSin {
public:
    enum SinSource {
        Prec,
        Quart,
        Full,
        Fixed,
        Lib,
    };
#if __cplusplus >= 202002L
    constexpr static const std::array<std::string, 5> sSinSource{"Prec", "Quart", "Full", "Fixed", "Lib"};
    constexpr static const std::string sinSourceName(SinSource sinSouce) { return sSinSource[static_cast<size_t>(sinSouce)]; }
#if __cplusplus >= 202302L
    constexpr
#endif
    static SinSource toSinSource(const std::string & name) {
        for (size_t i = 0; i < sSinSource.size(); ++i) {
            if (name == sSinSource[i])
                return static_cast<SinSource>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SinSource";
        throw std::invalid_argument(oss.str());
    }
#else
    static const char* const sSinSource[5];
    static std::string sinSourceName(SinSource sinSouce) { return sSinSource[static_cast<int>(sinSouce)]; }
    static SinSource toSinSource(const std::string & name) {
        for (size_t i = 0; i < 5; ++i) {
            if (name == sSinSource[i])
                return static_cast<SinSource>(i);
        }
        std::ostringstream oss;
        oss << "Invalid name=" << name << " for SinSource";
        throw std::invalid_argument(oss.str());
    }
#endif

#if __cplusplus >= 201103L
    constexpr static unsigned char MIN_ORDER = 1;
    constexpr static unsigned char MAX_ORDER = 18;
#else
    static const unsigned char MIN_ORDER;
    static const unsigned char MAX_ORDER;
#endif
        // (1 << MAX_ORDER) == PI
    static void validateOrder(unsigned char order);
        // assert MIN_ORDER <= order <= MAX_ORDER, and throw invalid_argument otherwise
    static const std::vector<VarDbl>& validateSinSource(SinSource sinSource);
        // get the sin source according to {sinSource}, and throw invalid_argument otherwise
    static unsigned char getOrder(size_t size);
        // get the order so that {size} = 2^{order}, and throw invalid_argument otherwise
    VarDbl sin(long long freq, unsigned char order) const;
        // sin(pi * freq / (1 << order))
    VarDbl cos(long long freq, unsigned char order) const;
        // cos(pi * freq / (1 << order))
    long long getIndex(long long freq, unsigned char order) const;

    const SinSource sinSource;

    bool dump(unsigned char order, const std::string& dumpPath) const;

    IndexSin(SinSource sinSource = Quart, const std::string& dumpDir = "");
        // if {dumpDir} is not empty, read value from the file in the {dumpDir}
private:
    const std::vector<VarDbl>& _sSin;

    static std::vector<VarDbl> _sSinQuart;
    static std::vector<VarDbl> _sSinFull;
    static std::vector<VarDbl> _sSinFixed;
    static std::vector<VarDbl> _sSinLib;
    static std::vector<VarDbl> _sSinPrec;

    static void read(std::vector<VarDbl>& sSin, const SinSource sinSource, const std::string& dumpDir);
        // read the sin values from the file in the {dumpDir}, and throw invalid_argument on error

};

std::vector<VarDbl> IndexSin::_sSinQuart;
std::vector<VarDbl> IndexSin::_sSinFull;
std::vector<VarDbl> IndexSin::_sSinFixed;
std::vector<VarDbl> IndexSin::_sSinLib;
std::vector<VarDbl> IndexSin::_sSinPrec;

#if __cplusplus < 202002L
const char* const IndexSin::sSinSource[5] = {"Prec", "Quart", "Full", "Fixed", "Lib"};
#endif
#if __cplusplus < 201103L
const unsigned char IndexSin::MIN_ORDER = 1;
const unsigned char IndexSin::MAX_ORDER = 18;
#endif

inline const std::vector<VarDbl>& IndexSin::validateSinSource(SinSource sinSource)
{
    switch (sinSource) {
    case Quart:
        return _sSinQuart;
    case Full:
        return _sSinFull;
    case Fixed:
        return _sSinFixed;
    case Lib:
        return _sSinLib;
    case Prec:
        return _sSinPrec;
    default:
        std::ostringstream oss;
        oss << "Unknown SinSource " <<  sinSource << " for fft.sin()";
        throw std::invalid_argument(oss.str());
    }
}


inline void IndexSin::validateOrder(unsigned char order)
{
    if (order < MIN_ORDER) {
        std::ostringstream oss;
        oss << "The order " << order << " < " << MIN_ORDER << " for fft.sin()";
        throw std::invalid_argument(oss.str());
    }
    if (MAX_ORDER < order) {
        std::ostringstream oss;
        oss << "The order " << order << " > " << MAX_ORDER << " for fft.sin()";
        throw std::invalid_argument(oss.str());
    }
}


inline unsigned char IndexSin::getOrder(size_t size)
{
    unsigned char order = IndexSin::MIN_ORDER;
    for (; order <= IndexSin::MAX_ORDER; ++order) {
        if ((1 << order) == (int)size) {
            return order;
        }
    }
    std::ostringstream oss;
    oss << "Invalid input array size " << size << " for fft.transform()";
    throw std::invalid_argument(oss.str());
}


inline long long IndexSin::getIndex(long long freq, unsigned char order) const
{
    validateOrder(order);
    bool pos = (freq >= 0);
    const size_t size = 1 << order;
    lldiv_t res = std::div(std::abs(freq), (long long)size);
    if (res.quot & 1)
        pos = !pos;
    if (((sinSource == Quart) || (sinSource == Prec)) && (res.rem > (long long)(size >> 1)))
        res.rem = size - res.rem;
    return pos? res.rem : - res.rem;
}


inline bool IndexSin::dump(unsigned char order, const std::string& dumpPath) const
{
    if (sinSource == Lib) {
        std::cerr << "Cannot dump sin values for SinSource::Lib";
        return false;
    }
    std::ofstream ofs(dumpPath.c_str());
    if (!ofs.is_open())
        return false;
    ofs << std::scientific << std::setprecision(20);
    ofs << "Index\tX\tValue\tUncertainty\n";
    const size_t size = 1 << order;
    for (size_t i = 0; i < _sSin.size(); ++i) {
        ofs << i << '\t' << ((double) i) / size << '\t' << _sSin[i].value() << '\t' << _sSin[i].uncertainty() << "\n";
   }
   return true;
}

inline void IndexSin::read(std::vector<VarDbl>& sSin, const SinSource sinSource, const std::string& dumpDir)
{
    if (sinSource == Lib) {
        throw std::invalid_argument("Cannot dump sin values for SinSource::Lib");
    }
    std::ostringstream oss;
    oss << dumpDir << "/IndexSin_" << sinSourceName(sinSource) << "_" << (unsigned) MAX_ORDER << ".txt";
    const std::string dumpPath = oss.str();
    oss.str("");
    std::ifstream ifs(dumpPath.c_str());
    if (!ifs.is_open()) {
        oss << "Failed to open file " << dumpPath;
        throw std::invalid_argument(oss.str());
    }
    std::string line;
    std::getline(ifs, line);
    if (line != "Index\tX\tValue\tUncertainty") {
        oss << "Invalid header in file " << dumpPath << ": " << line;
        throw std::invalid_argument(oss.str());
    }
    const size_t size = (1 << MAX_ORDER);
    sSin.reserve(size + 1);
    for (size_t n = 1; std::getline(ifs, line); ++n) {
        std::istringstream iss(line);
        int i;
        double x, value, uncertainty;
        if ((iss >> i >> x >> value >> uncertainty).fail()) {
            oss << "At line #" << n <<", failed to read value/uncertainty: " << line;
            throw std::invalid_argument(oss.str());
        }
#if __cplusplus >= 201103L
        sSin.emplace_back(value, uncertainty);
#else
        sSin.push_back(VarDbl(value, uncertainty));
#endif
    }
    if (sSin.size() != (size + 1)) {
        oss << "Only read " << sSin.size() << " lines from " << dumpPath << ", need " << (size + 1) << " lines";
        throw std::invalid_argument(oss.str());
    }
}




inline IndexSin::IndexSin(SinSource sinSource, const std::string& dumpDir) :
    sinSource(sinSource), _sSin(validateSinSource(sinSource))
{
    const size_t size = (1 << MAX_ORDER);
    const size_t half = size / 2;
    const size_t quart = size / 4;
    if (_sSinPrec.empty()) {
        if (dumpDir.empty()) {
            Histogram<double, size_t> histo(1.0, 20);
            for (size_t i = 0; i < size; ++i) {
                const double value = std::sin(_VAR_DBL_PI * i /size);
                if (i <= quart) {
                    const long double val = std::sin(_VAR_DBL_PI_LD * i /size);
#if __cplusplus >= 201103L
                    _sSinPrec.emplace_back((double)val, std::abs(val - ((double) val)));
                    _sSinQuart.emplace_back(value);
#else
                    _sSinPrec.push_back(VarDbl((double)val, std::abs(val - ((double) val))));
                    _sSinQuart.push_back(VarDbl(value));
#endif
                    const double ulpVal = ulp((double) val);
                    if (ulpVal > 0)
                        histo.addAt((val - ((double) val)) /ulpVal, i);
                }  else if (i <= half) {
                    const double value2 = std::cos(_VAR_DBL_PI * (half - i) /size);
                    const long double val = std::cos(_VAR_DBL_PI_LD * (half - i) /size);
#if __cplusplus >= 201103L
                    _sSinPrec.emplace_back((double)val, std::abs(val - value2));
                    _sSinQuart.emplace_back(value2);
#else
                    _sSinPrec.push_back(VarDbl((double)val, std::abs(val - value2)));
                    _sSinQuart.push_back(VarDbl(value2));
#endif
                    const double ulpVal = ulp((double) val);
                    if (ulpVal > 0)
                        histo.addAt((val - ((double) val)) /ulpVal, half - i);
                }
#if __cplusplus >= 201103L
                _sSinFull.emplace_back(value);
                _sSinFixed.emplace_back(value, VarDbl::ulp(1.));
#else
                _sSinFull.push_back(VarDbl(value));
                _sSinFixed.push_back(VarDbl(value, VarDbl::ulp(1.)));
#endif
            }
            assert(_sSinPrec.size() == (half + 1));
            assert(_sSinQuart.size() == (half + 1));
            assert(_sSinFull.size() == size);
            assert(_sSinFixed.size() == size);
            std::ofstream ofs("Output/Prec.histo.txt");
            ofs << "Count\tMean\tDev\tMin\tMinAt\tMax\tMaxAt\tLower Count\tUpper Count" << histo.header() << "\n";
            ofs << histo.count() << '\t' << histo.mean() << '\t' << histo.std()
                << '\t' << histo.min() << '\t' << histo.minAt().value()
                << '\t' << histo.max() << '\t' << histo.maxAt().value()
                << '\t' << histo.lowers() << '\t' << histo.uppers() << histo.formatted();
            ofs.flush();
        } else {
            read(_sSinPrec, Prec, dumpDir);
            read(_sSinQuart, Quart, dumpDir);
            read(_sSinFull, Full, dumpDir);
            read(_sSinFixed, Fixed, dumpDir);
            _sSinPrec.erase(_sSinQuart.begin() + half + 1, _sSinQuart.end());
            _sSinQuart.erase(_sSinPrec.begin() + half + 1, _sSinPrec.end());
            _sSinFull.erase(_sSinFull.begin() + size, _sSinFull.end());
            _sSinFixed.erase(_sSinFixed.begin() + size, _sSinFixed.end());
        }
        assert(_sSinLib.empty());
    }
}

inline VarDbl IndexSin::sin(long long freq, unsigned char order) const
{
    validateOrder(order);
    if (sinSource == Lib) {
        const double value = std::sin(_VAR_DBL_PI * freq / (1 << order));
        return VarDbl(value);
    }
    const long long idx = getIndex(freq, order);
    const VarDbl val = _sSin[abs(idx) << (MAX_ORDER - order)];
    return (0 <= idx)? val : -val;
}

inline VarDbl IndexSin::cos(long long freq, unsigned char order) const
{
    return sin(freq + (1 << (order - 1)), order);
}




} // namespace var_dbl

#undef _VAR_DBL_PI
#undef _VAR_DBL_PI_LD

#endif  // __IndexSin_h__
