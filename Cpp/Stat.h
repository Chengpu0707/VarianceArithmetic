#include <algorithm>
#include <cmath>
#include <sstream>
#include <random>
#include <format>


#ifndef __Stat_h__
#define __Stat_h__
namespace var_dbl 
{

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
        std::normal_distribution d{mean, dev};
        return d(gen);
    }

    double white() {
        std::uniform_real_distribution d{mean - std::sqrt(3)* dev, mean + std::sqrt(3)* dev};
        return d(gen);
    }
};


template <typename T, typename HINT = int> requires std::floating_point<T> || std::integral<T>
class Stat {
    T _min;
    T _max;
    HINT _minAt;
    HINT _maxAt;
    bool _at;
    unsigned _count;
    double _sum;
    double _sum2;

public: 
    Stat() { clear(); }
    template<typename U> Stat(std::initializer_list<U> sInit);
    template<typename InputIt> Stat(InputIt begin, InputIt end);

    void clear();

    unsigned add(const T val);
    unsigned addAt(const T val, HINT at);
    template<typename U> unsigned add(std::initializer_list<U> sInit);
    template<typename InputIt> unsigned add(InputIt begin, InputIt end);

    unsigned count() const { return _count; }
    T min() const { return _min; }
    T max() const { return _max; }
    HINT minAt() const;
    HINT maxAt() const;
    double mean() const;
    double std() const;
};


class Histogram : public Stat<double, int> 
{
	const unsigned _center;
	std::vector<unsigned> _sHistogram;
    unsigned _lowers = 0, _uppers = 0;

public:
	const double range;
	const unsigned divids;
    unsigned lowers() const { return _lowers; }
    unsigned uppers() const { return _uppers; }

    Histogram(double range = 3, unsigned divids = 5) :
        Stat<double>(), range(range), divids(divids), _center(divids * range)
    {
        const unsigned size = 1 + (_center * 2);
        _sHistogram.reserve(size);
        _sHistogram.insert(_sHistogram.begin(), size, 0);
    }

    void clear();
    unsigned add(double val);
    template<typename U> unsigned add(std::initializer_list<U> sInit);
    template<typename InputIt> unsigned add(InputIt begin, InputIt end);

    std::vector<unsigned> histogram() { return std::vector<unsigned>(_sHistogram); }
    std::string header();
    std::string formatted(bool normalized = true);
};


template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline HINT Stat<T, HINT>::minAt() const 
{ 
    if (_at) 
        return _minAt; 
    throw std::domain_error("minAt(): No at info"); 
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline HINT Stat<T, HINT>::maxAt() const 
{ 
    if (_at) 
        return _maxAt; 
    throw std::domain_error("maxAt(): No at info"); 
}


template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline void Stat<T, HINT>::clear() 
{
    _min = std::numeric_limits<T>::max();
    _max = -std::numeric_limits<T>::max();
    _at = false;
    _count = 0;
    _sum = 0;
    _sum2 = 0;
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
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

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline unsigned Stat<T, HINT>::addAt(const T val, HINT at) 
{
    _at = true;
    if (_min > val)
        _minAt = at;
    if (_max < val)
        _maxAt = at;
    return add(val);
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
template<typename U> 
inline unsigned Stat<T, HINT>::add(std::initializer_list<U> sInit)
{
    for (auto v: sInit)
        add(v);
    return _count;
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
template <typename InputIt>
inline unsigned Stat<T, HINT>::add(InputIt begin, InputIt end) 
{
    for (auto it = begin; it != end; ++it)
        add(*it);
    return _count;
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
template<typename InputIt> 
inline Stat<T, HINT>::Stat(InputIt begin, InputIt end) 
{
    clear();
    add(begin, end);
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
template<typename U> 
inline Stat<T, HINT>::Stat(std::initializer_list<U> sInit) 
{
    clear();
    for (auto v: sInit)
        add(v);
}


template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline double Stat<T, HINT>::mean() const 
{
    if (_count <= 0)
        return std::numeric_limits<double>::quiet_NaN();
    return _sum / _count;
}


template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline double Stat<T, HINT>::std() const 
{
    if (_count <= 0)
        return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(std::max((double) 0, (_sum2 - (_sum*_sum / _count)) / (_count - 1)));
}


inline void Histogram::clear()
{
    _lowers = 0;
    _uppers = 0;
    std::fill(_sHistogram.begin(), _sHistogram.end(), 0);
    Stat<double>::clear();
}


inline unsigned Histogram::add(double val)
{
    const int idx = (int) std::round( val * divids ) + _center;
    if (idx < 0)
        ++_lowers;
    else if (idx >= _sHistogram.size())
        ++_uppers;
    else
        ++_sHistogram[idx];
    return Stat<double>::add(val);
}


template<typename U> inline unsigned Histogram::add(std::initializer_list<U> sInit)
{
    for (auto v: sInit)
        add(v);
    return count();
}


template<typename InputIt> inline unsigned Histogram::add(InputIt begin, InputIt end)
{
    for (auto it = begin; it != end; ++it)
        add(*it);
    return count();
}


inline std::string Histogram::header()
{
    std::ostringstream oss;
    for (double i = -(int) _center; i <= (int) _center; ++i)
        oss << "\t" << i / divids;
    return oss.str();
}


inline std::string Histogram::formatted(bool normalized)
{
    const double cnt = count() - lowers() - uppers();
    std::ostringstream oss;
    for (size_t i = 0; i < _sHistogram.size(); ++i)
        oss << "\t" << (((cnt == 0) || !normalized)? _sHistogram[i] : _sHistogram[i]/cnt);
    return oss.str();
}


}   // var_dbl
#endif  // __Stat_h__