/*
Statistical functions
*/

#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

#ifndef __Stat_h__
#define __Stat_h__
namespace var_dbl 
{

struct Stat {
    double mean;
    double stddev;
    unsigned count;
};

template< class InputIt >
inline Stat calcStat(InputIt begin, InputIt end) 
{
    Stat res;
    res.count = end - begin;
    res.mean = std::nan("1");
    res.stddev = std::nan("2");
    if (res.count <= 0)
        return res;
    res.mean = std::accumulate(begin, end, 0.0) / res.count;
    if (res.count <= 1)
        return res;
    res.stddev = std::accumulate(begin, end, 0, [res](double init, double item){
            return init + (item - res.mean) * (item - res.mean);
        }) / (res.count - 1);
    res.stddev = std::sqrt(res.stddev);
    return res;
}

struct Histo {
    Stat stat;
    double range;
    unsigned divides;
    unsigned less;
    unsigned more;
    std::vector<std::pair<float, float>> sHisto;
};

template< class InputIt >
inline Histo calcHisto(InputIt begin, InputIt end, double range, unsigned divides,
                        bool normalized=true) {
    Histo res;
    res.range = range;
    res.divides = divides;
    res.less = 0;
    res.more = 0;
    res.stat = calcStat(begin, end);
    if (!std::isfinite(res.stat.stddev) || (range <= 0) || (divides == 0))
        return res;
    int size = std::round(range * divides);
    res.range = (double) size / divides;
    for (int i = -size; i < size; ++i) 
        res.sHisto.push_back({(i + 0.5f)/divides, 0});
    for (auto it = begin; it != end; ++it) {
        int idx = std::floor((*it - res.stat.mean) / res.stat.stddev * divides);
        if (idx < -size)
            ++res.less;
        else if (idx >= size)
            ++res.more;
        else 
            ++res.sHisto[idx + size].second;
    }
    const int count = res.stat.count - res.less - res.more;
    if (normalized && (count > 0)) {
        for (int i = -size; i <= size; ++i) 
            res.sHisto[i + size].second /= count;
    }
    return res;
}


} // namespace var_dbl
#endif  // __Stat_h__