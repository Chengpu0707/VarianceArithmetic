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

class Momentum {
    const std::vector<double> _sMomentum;

protected:
    Momentum(const std::vector<double>& sMomentum, double bounding, double leakage) :
        _sMomentum(sMomentum), bounding(bounding), leakage(leakage) {
    }

public:
    const double bounding;
    const double leakage;
    double operator[](size_t i) const { return _sMomentum[i]; }
    size_t maxOrder() const { return _sMomentum.size(); }
};


class NormalMomentum : public Momentum {
    std::vector<double> init(double bounding);

public:
    static double pdf(double z);

    NormalMomentum(double bounding);
};


inline double NormalMomentum::pdf(double z) {
    const double x2 = z * z;
    return std::exp(-0.5 * x2) / sqrt(2 * std::numbers::pi);
}



inline NormalMomentum::NormalMomentum(double bounding) 
        : Momentum(NormalMomentum::init(bounding), bounding, 1 - std::erf(bounding/std::sqrt(2))) { 
}

inline std::vector<double> NormalMomentum::init(double bounding) {
    const size_t maxOrder = 10000;
    const double bounding2 = bounding * bounding;
    double term = 2 * pdf(bounding) * bounding;
    double sTerm[maxOrder];
    size_t n = 0;
    for (; n < maxOrder; ++n) {
        sTerm[n] = term / (2*n + 1);
        if (!std::isfinite(sTerm[n])) {
            break;
        }
        term *= bounding2;
    }

    std::vector<double> sMomentum(n << 1, 0);
    for (size_t i = 0; i < n; ++i) 
        sMomentum[i << 1] = sTerm[i];
    for (size_t j = 2; j < sMomentum.size(); ++j) {
        for (size_t i = 0; i < n; ++i) {
            sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2;
            const double prev = sMomentum[i << 1];
            sMomentum[i << 1] += sTerm[i];
            if (prev == sMomentum[i << 1]) {
                n = i;
                break;
            }
        }
        if (n <= 0) 
            break;
    }
    return sMomentum;
}



} // namespace var_dbl
#endif // __Momentum_h__