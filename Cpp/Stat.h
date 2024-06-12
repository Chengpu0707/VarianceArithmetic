/*
Statistical functions
*/

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

#ifndef __Stat_h__
#define __Stat_h__
namespace var_dbl 
{

struct Stat {
    double min = std::numeric_limits<double>::max();
    double max = -std::numeric_limits<double>::max();
    double stddev = 0;
    double mean = 0;
    unsigned count = 0;
};

template< class InputIt >
Stat calcStat(InputIt begin, InputIt end);

struct Histo {
    Stat stat;
    float range = 0;       // centered at mean and normalized by stddv, may be different from the input
    unsigned divides = 0;  // per stddev
    unsigned less = 0;     // less than range count
    unsigned more = 0;      // more than range count
    std::vector<float> sHisto;  // 2*range*divides buckets around stat.mean
                                // center for each bucket is by calcHistoCenters()
};

std::vector<float> calcHistoCenters(float range, unsigned divides);

template< class InputIt >
bool calcHisto(Histo& res, InputIt begin, InputIt end, float range, unsigned divides, bool normalized=true);




template< class InputIt >
inline Stat calcStat(InputIt begin, InputIt end) 
{
    Stat res;
    res.count = end - begin;
    if (res.count <= 0)
        return res;
    res.mean = std::accumulate(begin, end, 0.0) / res.count;
    for (auto it = begin; it != end; ++it) {
        if (res.min > *it)
            res.min = *it;
        if (res.max < *it)
            res.max = *it;
        res.stddev += (*it - res.mean) * (*it - res.mean);
    }
    if (res.count <= 1) {
        return res;
    }
    res.stddev = std::sqrt(res.stddev/(res.count - 1));
    return res;
}


inline std::vector<float> calcHistoCenters(float range, unsigned divides)
{
    std::vector<float> sRes;
    if ((range <= 0) || (divides == 0))
        return sRes;
    int size = std::round(range * divides);
    for (int i = -size; i < size; ++i) 
        sRes.push_back((i + 0.5f)/divides);
    return sRes;
}

template< class InputIt >
inline bool calcHisto(Histo& res, InputIt begin, InputIt end, float range, unsigned divides, bool normalized) {
    res.range = range;
    res.divides = divides;
    res.less = 0;
    res.more = 0;
    res.stat = calcStat(begin, end);
    if (!std::isfinite(res.stat.stddev) || (range <= 0) || (divides == 0))
        return false;
    int size = std::round(range * divides);
    res.range = (float) size / divides;
    res.sHisto.clear();
    res.sHisto.reserve(2*size + 1);
    res.sHisto.insert(res.sHisto.end(), 2*size, 0);
    for (auto it = begin; it != end; ++it) {
        int idx = std::floor((*it - res.stat.mean) / res.stat.stddev * divides);
        if (idx < -size)
            ++res.less;
        else if (idx >= size)
            ++res.more;
        else 
            ++res.sHisto[idx + size];
    }
    const int count = res.stat.count - res.less - res.more;
    if (normalized && (count > 0)) {
        for (int i = -size; i <= size; ++i) 
            res.sHisto[i + size] /= count;
    }
    return true;
}


} // namespace var_dbl
#endif  // __Stat_h__