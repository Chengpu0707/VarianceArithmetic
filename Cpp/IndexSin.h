/*
A sin() with index frequence as input, with resolution \pi/size()
size() = 1 << "order".
When the index "freq" == size, the result is sin(\pi).
This function has minimal float error when "withUncertainty"==False.

When "withUncertainty"==True, use regression to calculate the sin

*/

#include "VarDbl.h"

#include <cmath>
#include <sstream>
#include <vector>


#ifndef __IndexSin__
#define __IndexSin__
namespace var_dbl 
{

class IndexSin {
    const size_t _size;
    const int _half;
    std::vector<VarDbl> _sSin;

public:
    constexpr const static double PI = 3.14159265358979323846;

    IndexSin(unsigned order, bool withUncertainty);
    size_t size() const { return _size; }

    VarDbl sin(int freq) const;
    VarDbl cos(int freq) const;
    VarDbl tan(int freq) const;

    int get_sin_index(int freq) const;
        // get index into _sSin, with -index means -sin

};


inline IndexSin::IndexSin(unsigned order, bool withUncertainty) :
    _size(1 << order), _half(1 << (order - 1))
{
    if (order < 3) {
        std::ostringstream os;
        os << "order=" << order << " < 3 for IndexSin.";
        throw std::invalid_argument(os.str());
    }
    if (withUncertainty) {
        std::ostringstream os;
        os << "withUncertainty=" << withUncertainty<< " not implemented yet for IndexSin.";
        throw std::invalid_argument(os.str());
    }
    _sSin.reserve(_size + 1);
    if (! withUncertainty) {
        for (int i = 0; i < _half/2; ++i)
            _sSin.push_back(std::sin(PI/_size*i));
        for (int i = 0; i <= _half/2; ++i)
            _sSin.push_back(std::cos(PI*1/4 - PI*i/_size));
        return;
    }
}

inline int IndexSin::get_sin_index(int freq) const
{
    int div = freq / _half, rem = freq % _half;
    if (div & 1) {
        if (div > 0) {
            div -= 1;
            rem = _half - rem;
        } else {
            div += 1;
            rem = -_half - rem;
        }
    }
    if (div & 2)
        rem = -rem;
    return rem;
}


inline VarDbl IndexSin::sin(int freq) const
{
    const int idx = get_sin_index(freq);
    return (idx >= 0)? _sSin[idx] : -_sSin[-idx];
}

inline VarDbl IndexSin::cos(int freq) const
{
    return sin(freq + _half);
}

inline VarDbl IndexSin::tan(int freq) const
{
//TODO:    return sin(freq) / cos(freq);
}



} // namespace testvar_dbl
#endif  // __IndexSin__