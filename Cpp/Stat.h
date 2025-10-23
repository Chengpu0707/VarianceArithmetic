#include <algorithm>
#include <cmath>
#include <optional>
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
protected:
    T _min;
    T _max;
    std::optional<HINT> _minAt;
    std::optional<HINT> _maxAt;
    unsigned _count;
    double _sum;
    double _sum2;

public: 
    Stat() { clear(); }
    Stat(std::initializer_list<T> sInit);
    template<typename InputIt> Stat(InputIt begin, InputIt end);

    void clear();

    unsigned add(const T val);
    unsigned addAt(const T val, HINT at);
    unsigned add(std::initializer_list<T> sInit);
    template<typename InputIt> unsigned add(InputIt begin, InputIt end);

    unsigned count() const { return _count; }
    T min() const { return _min; }
    T max() const { return _max; }
    std::optional<HINT> minAt() const { return _minAt; }
    std::optional<HINT> maxAt() const { return _maxAt; }
    double mean() const;
    double std() const;
};

template <typename T = double, typename HINT = int> requires std::floating_point<T>
class Histogram : public Stat<T, HINT> 
{
protected:
	const unsigned _center;
	std::vector<unsigned> _sHistogram;
    unsigned _lowers = 0, _uppers = 0;

public:
	const T range;
	const unsigned divids;
    unsigned lowers() const { return _lowers; }
    unsigned uppers() const { return _uppers; }

    Histogram(T range = 3, unsigned divids = 5);

    void clear();
    unsigned add(T val);
    unsigned addAt(T val, HINT at);
    unsigned add(std::initializer_list<T> sInit);
    template<typename InputIt> unsigned add(InputIt begin, InputIt end);

    std::vector<unsigned> histogram() { return std::vector<unsigned>(_sHistogram); }
    std::string header() const;
    std::string formatted(bool normalized = true) const;
};


template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
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
    if (_min > val)
        _minAt = at;
    if (_max < val)
        _maxAt = at;
    return add(val);
}

template <typename T, typename HINT> requires std::floating_point<T> || std::integral<T>
inline unsigned Stat<T, HINT>::add(std::initializer_list<T> sInit)
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
inline Stat<T, HINT>::Stat(std::initializer_list<T> sInit) 
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


template <typename T, typename HINT> requires std::floating_point<T>
inline Histogram<T, HINT>::Histogram(T range, unsigned divids) :
        Stat<T, HINT>(), range(range), divids(divids), _center(divids * range)
{
    const unsigned size = 1 + (_center * 2);
    _sHistogram.reserve(size);
    _sHistogram.insert(_sHistogram.begin(), size, 0);
}


template <typename T, typename HINT> requires std::floating_point<T>
inline void Histogram<T, HINT>::clear()
{
    _lowers = 0;
    _uppers = 0;
    std::fill(_sHistogram.begin(), _sHistogram.end(), 0);
    Stat<T, HINT>::clear();
}


template <typename T, typename HINT> requires std::floating_point<T>
inline unsigned Histogram<T, HINT>::add(T val)
{
    const int idx = (int) std::round( val * divids ) + _center;
    if (idx < 0)
        ++_lowers;
    else if (idx >= _sHistogram.size())
        ++_uppers;
    else
        ++_sHistogram[idx];
    return Stat<T, HINT>::add(val);
}


template <typename T, typename HINT> requires std::floating_point<T>
inline unsigned Histogram<T, HINT>::addAt(T val, HINT at)
{
    if (Stat<T, HINT>::_min > val)
        Stat<T, HINT>::_minAt = at;
    if (Stat<T, HINT>::_max < val)
        Stat<T, HINT>::_maxAt = at;
    return Histogram<T, HINT>::add(val);
}


template <typename T, typename HINT> requires std::floating_point<T>
inline unsigned Histogram<T, HINT>::add(std::initializer_list<T> sInit)
{
    for (auto v: sInit)
        add(v);
    return Stat<T, HINT>::count();
}


template <typename T, typename HINT> requires std::floating_point<T>
template<typename InputIt> 
inline unsigned Histogram<T, HINT>::add(InputIt begin, InputIt end)
{
    for (auto it = begin; it != end; ++it)
        add(*it);
    return Stat<T, HINT>::count();
}


template <typename T, typename HINT> requires std::floating_point<T>
inline std::string Histogram<T, HINT>::header() const
{
    std::ostringstream oss;
    for (double i = -(int) _center; i <= (int) _center; ++i)
        oss << "\t" << i / divids;
    return oss.str();
}


template <typename T, typename HINT> requires std::floating_point<T>
inline std::string Histogram<T, HINT>::formatted(bool normalized) const
{
    const double cnt = Stat<T, HINT>::count() - lowers() - uppers();
    std::ostringstream oss;
    for (size_t i = 0; i < _sHistogram.size(); ++i)
        oss << "\t" << (((cnt == 0) || !normalized)? _sHistogram[i] : _sHistogram[i]/cnt);
    return oss.str();
}


}   // var_dbl
#endif  // __Stat_h__