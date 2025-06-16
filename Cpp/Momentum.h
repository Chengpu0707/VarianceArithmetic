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
    static double pdf(double z);

    NormalMomentum(double BOUNDING=5, size_t MAX_ORDER = 10000);
    double operator[](size_t i) const;
    size_t maxOrder() const { return _sMomentum.size() * 2; }

    const double BOUNDING;
    const double LEAKAGE;
private:
    std::vector<double> _sMomentum;
};


inline double NormalMomentum::pdf(double z) {
    const double x2 = z * z;
    return std::exp(-0.5 * x2) / sqrt(2 * std::numbers::pi);
}



inline NormalMomentum::NormalMomentum(double BOUNDING, size_t MAX_ORDER) 
        : BOUNDING(BOUNDING), LEAKAGE(1 - std::erf(BOUNDING/std::sqrt(2))) { 

    double term = 2 * pdf(BOUNDING) * BOUNDING;
    const double bounding2 = BOUNDING * BOUNDING;
    double sTerm[MAX_ORDER];
    size_t n = 0;
    for (; n < MAX_ORDER; ++n) {
        sTerm[n] = term / (2*n + 1);
        if (!std::isfinite(sTerm[n])) {
            break;
        }
        term *= bounding2;
    }
    
    _sMomentum.insert(_sMomentum.end(), sTerm, sTerm + n);
    for (size_t j = 2; j < _sMomentum.size(); ++j) {
        for (size_t i = 0; i < n; ++i) {
            sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2;
            const double prev = _sMomentum[i];
            _sMomentum[i] += sTerm[i];
            if (prev == _sMomentum[i]) {
                n = i;
                break;
            }
        }
        if (n <= 0) 
            break;
    }
}


inline double NormalMomentum::operator[](size_t i) const {
    const size_t j = i >> 1;
    if (j >= _sMomentum.size() || (i & 1) != 0)
        return 0;
    return _sMomentum[j];
}

} // namespace var_dbl
#endif // __Momentum_h__