/*
NormalMomentum is needed for Taylor expansion in VarDbl.  
It is a constexpr class.
*/

#include <cmath>
#include <numbers>
#include <vector>


#ifndef __Momentum_h__
#define __Momentum_h__
namespace var_dbl 
{

class NormalMomentum {
public:
    NormalMomentum(double BINDING=5, size_t DIVID=128);
    double operator[](size_t i) const;
    size_t size() const { return _sFactor.size(); }

    const double BINDING;
private:
    std::vector<double> _sFactor;
};


inline NormalMomentum::NormalMomentum(double BINDING, size_t DIVID) 
        : BINDING(BINDING), _sFactor(2000) { 
    const double divid2 = DIVID * DIVID;
    const double norm = 1.0/sqrt(2*std::numbers::pi) / DIVID;
    const long limit = DIVID * BINDING;
    for (long i = -limit ; i < limit; ++i) {
        const double x2 = (i + 0.5)*(i + 0.5) / divid2;
        const double pdf = norm * exp(- x2 * 0.5);
        double sq = 1;
        for (int j = 0; j < _sFactor.size(); j += 2) {
            if (! std::isfinite(_sFactor[j]))
                break;
            _sFactor[j] += pdf * sq;
            sq *= x2;
        }
    }
    size_t i;
    for (i = 0; i < _sFactor.size(); ++i) {
        if (! std::isfinite(_sFactor[i]))
            break;
    }
    _sFactor.erase(_sFactor.begin() + i, _sFactor.end());
}


inline double NormalMomentum::operator[](size_t i) const {
    if (i >= _sFactor.size())
        return 0;
    return _sFactor[i];
}

} // namespace var_dbl
#endif // __Momentum_h__