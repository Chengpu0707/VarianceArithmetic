#include <algorithm>
#include <cmath>
#include <limits>
#if __cplusplus >= 202002L
#include <format>
#endif
#if __cplusplus >= 201703L
#include <optional>
#endif
#if __cplusplus >= 201103L
#include <random>
#else
#include <cstdlib>
#endif
#include <sstream>
#include <string>
#include <vector>


#ifndef __Stat_h__
#define __Stat_h__
#if __cplusplus >= 202002L
#define STAT_TMPL template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
#define HISTO_TMPL template <typename T, typename HINT> requires std::floating_point<T>
#else
#define STAT_TMPL template <typename T, typename HINT>
#define HISTO_TMPL template <typename T, typename HINT>
#endif
namespace var_dbl
{

#if __cplusplus < 201703L
template <typename T>
struct _Optional {
    bool _has;
    T _val;
    _Optional() : _has(false), _val() {}
    _Optional(const T& v) : _has(true), _val(v) {}
    bool has_value() const { return _has; }
    T value() const { return _val; }
    void reset() { _has = false; }
    _Optional& operator=(const T& v) { _has = true; _val = v; return *this; }
};
#endif


#if __cplusplus >= 201103L
class Random
{
    std::random_device rd{};
    std::mt19937 gen{rd()};

public:
    const double mean;
    const double dev;

    Random(double mean, double dev) :
        mean(mean), dev(std::abs(dev))
    {}

    double gauss() {
        if (dev == 0) return mean;
        std::normal_distribution<double> d(mean, dev);
        return d(gen);
    }

    double white() {
        if (dev == 0) return mean;
        std::uniform_real_distribution<double> d(mean - std::sqrt(3.0)* dev, mean + std::sqrt(3.0)* dev);
        return d(gen);
    }
};
#else
class Random
{
public:
    const double mean;
    const double dev;

    Random(double mean, double dev) :
        mean(mean), dev(std::abs(dev))
    {}

    double gauss() {
        // Box-Muller transform
        const double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
        const double u2 = rand() / (RAND_MAX + 1.0);
        const double pi = 3.14159265358979323846;
        return mean + dev * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * pi * u2);
    }

    double white() {
        const double lo = mean - std::sqrt(3.0) * dev;
        const double hi = mean + std::sqrt(3.0) * dev;
        return lo + (hi - lo) * (rand() / (RAND_MAX + 1.0));
    }
};
#endif


#if __cplusplus >= 202002L
template <typename T, typename HINT = int> requires std::floating_point<T> || std::integral<T>
#else
template <typename T, typename HINT = int>
#endif
class Stat {
protected:
    T _min;
    T _max;
#if __cplusplus >= 201703L
    std::optional<HINT> _minAt;
    std::optional<HINT> _maxAt;
#else
    _Optional<HINT> _minAt;
    _Optional<HINT> _maxAt;
#endif
    unsigned _count;
    double _sum;
    double _sum2;

public:
    Stat() { clear(); }
#if __cplusplus >= 201103L
    Stat(std::initializer_list<T> sInit);
#endif
    template<typename InputIt> Stat(InputIt begin, InputIt end);

    void clear();

    unsigned add(const T val);
    unsigned addAt(const T val, HINT at);
#if __cplusplus >= 201103L
    unsigned add(std::initializer_list<T> sInit);
#endif
    template<typename InputIt> unsigned add(InputIt begin, InputIt end);

    unsigned count() const { return _count; }
    T min() const { return _min; }
    T max() const { return _max; }
#if __cplusplus >= 201703L
    std::optional<HINT> minAt() const { return _minAt; }
    std::optional<HINT> maxAt() const { return _maxAt; }
#else
    _Optional<HINT> minAt() const { return _minAt; }
    _Optional<HINT> maxAt() const { return _maxAt; }
#endif
    double mean() const;
    double std() const;
};

#if __cplusplus >= 202002L
template <typename T = double, typename HINT = int> requires std::floating_point<T>
#else
template <typename T = double, typename HINT = int>
#endif
class Histogram : public Stat<T, HINT>
{
protected:
	const unsigned _center;
	std::vector<unsigned> _sHistogram;
#if __cplusplus >= 201103L
    unsigned _lowers = 0, _uppers = 0;
#else
    unsigned _lowers, _uppers;
#endif

public:
	const T range;
	const unsigned divids;
    unsigned lowers() const { return _lowers; }
    unsigned uppers() const { return _uppers; }

    Histogram(T range = 3, unsigned divids = 5);

    void clear();
    unsigned add(T val);
    unsigned addAt(T val, HINT at);
#if __cplusplus >= 201103L
    unsigned add(std::initializer_list<T> sInit);
#endif
    template<typename InputIt> unsigned add(InputIt begin, InputIt end);

    std::vector<unsigned> histogram() { return std::vector<unsigned>(_sHistogram); }
    std::string header() const;
    std::string formatted(bool normalized = true) const;
};


STAT_TMPL
inline void Stat<T, HINT>::clear()
{
    _min = std::numeric_limits<T>::max();
    _max = -std::numeric_limits<T>::max();
    _minAt.reset();
    _maxAt.reset();
    _count = 0;
    _sum = 0;
    _sum2 = 0;
}

STAT_TMPL
inline unsigned Stat<T, HINT>::add(const T val)
{
    ++_count;
    if (_min > val)
        _min = val;
    if (_max < val)
        _max = val;
    _sum += val;
    _sum2 += val * val;
    return _count;
}

STAT_TMPL
inline unsigned Stat<T, HINT>::addAt(const T val, HINT at)
{
    if (_min > val)
        _minAt = at;
    if (_max < val)
        _maxAt = at;
    return add(val);
}

#if __cplusplus >= 201103L
STAT_TMPL
inline unsigned Stat<T, HINT>::add(std::initializer_list<T> sInit)
{
    for (auto v: sInit)
        add(v);
    return _count;
}
#endif

STAT_TMPL
template <typename InputIt>
inline unsigned Stat<T, HINT>::add(InputIt begin, InputIt end)
{
#if __cplusplus >= 201103L
    for (auto it = begin; it != end; ++it)
#else
    for (InputIt it = begin; it != end; ++it)
#endif
        add(*it);
    return _count;
}

STAT_TMPL
template<typename InputIt>
inline Stat<T, HINT>::Stat(InputIt begin, InputIt end)
{
    clear();
    add(begin, end);
}

#if __cplusplus >= 201103L
STAT_TMPL
inline Stat<T, HINT>::Stat(std::initializer_list<T> sInit)
{
    clear();
    for (auto v: sInit)
        add(v);
}
#endif


STAT_TMPL
inline double Stat<T, HINT>::mean() const
{
    if (_count <= 0)
        return std::numeric_limits<double>::quiet_NaN();
    return _sum / _count;
}


STAT_TMPL
inline double Stat<T, HINT>::std() const
{
    if (_count <= 0)
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(std::max((double) 0, (_sum2 - (_sum*_sum / _count)) / (_count - 1)));
}


HISTO_TMPL
inline Histogram<T, HINT>::Histogram(T range, unsigned divids) :
        Stat<T, HINT>(), range(range), divids(divids), _center(divids * range)
#if __cplusplus < 201103L
        , _lowers(0), _uppers(0)
#endif
{
    const unsigned size = 1 + (_center * 2);
    _sHistogram.reserve(size);
    _sHistogram.insert(_sHistogram.begin(), size, 0);
}


HISTO_TMPL
inline void Histogram<T, HINT>::clear()
{
    _lowers = 0;
    _uppers = 0;
    std::fill(_sHistogram.begin(), _sHistogram.end(), 0);
    Stat<T, HINT>::clear();
}


HISTO_TMPL
inline unsigned Histogram<T, HINT>::add(T val)
{
#if __cplusplus >= 201103L
    const int idx = (int) std::round( val * divids ) + _center;
#else
    const int idx = (int) ::round( val * divids ) + _center;
#endif
    if (idx < 0)
        ++_lowers;
    else if (idx >= (int)_sHistogram.size())
        ++_uppers;
    else
        ++_sHistogram[idx];
    return Stat<T, HINT>::add(val);
}


HISTO_TMPL
inline unsigned Histogram<T, HINT>::addAt(T val, HINT at)
{
    if (Stat<T, HINT>::_min > val)
        Stat<T, HINT>::_minAt = at;
    if (Stat<T, HINT>::_max < val)
        Stat<T, HINT>::_maxAt = at;
    return Histogram<T, HINT>::add(val);
}


#if __cplusplus >= 201103L
HISTO_TMPL
inline unsigned Histogram<T, HINT>::add(std::initializer_list<T> sInit)
{
    for (auto v: sInit)
        add(v);
    return Stat<T, HINT>::count();
}
#endif


HISTO_TMPL
template<typename InputIt>
inline unsigned Histogram<T, HINT>::add(InputIt begin, InputIt end)
{
#if __cplusplus >= 201103L
    for (auto it = begin; it != end; ++it)
#else
    for (InputIt it = begin; it != end; ++it)
#endif
        add(*it);
    return Stat<T, HINT>::count();
}


HISTO_TMPL
inline std::string Histogram<T, HINT>::header() const
{
    std::string result;
    for (double i = -(int) _center; i <= (int) _center; ++i) {
#if __cplusplus >= 202002L
        result += std::format("\t{}", i / divids);
#else
        std::ostringstream _oss; _oss << "\t" << (i / divids); result += _oss.str();
#endif
    }
    return result;
}


HISTO_TMPL
inline std::string Histogram<T, HINT>::formatted(bool normalized) const
{
    const double cnt = Stat<T, HINT>::count() - lowers() - uppers();
    std::string result;
    for (size_t i = 0; i < _sHistogram.size(); ++i) {
#if __cplusplus >= 202002L
        result += std::format("\t{}", ((cnt == 0) || !normalized) ? _sHistogram[i] : _sHistogram[i]/cnt);
#else
        std::ostringstream _oss;
        _oss << "\t" << (((cnt == 0) || !normalized) ? _sHistogram[i] : _sHistogram[i]/cnt);
        result += _oss.str();
#endif
    }
    return result;
}


}   // var_dbl
#endif  // __Stat_h__
