/*
NormalMoment is needed for Taylor expansion in VarDbl.
It is a constexpr class.
*/

#include <cmath>
#if __cplusplus >= 202002L
#include <numbers>
#endif
#include <vector>


#ifndef __Moment_h__
#define __Moment_h__
namespace var_dbl
{

class Moment {
    const std::vector<double> _sMoment;

protected:
    Moment(const std::vector<double>& sMoment, double bounding, double leakage) :
        _sMoment(sMoment), bounding(bounding), leakage(leakage) {
    }

public:
    const double bounding;
    const double leakage;
    double operator[](size_t i) const { return _sMoment[i]; }
    size_t maxOrder() const { return _sMoment.size(); }
};


class NormalMoment : public Moment {
    std::vector<double> init(double bounding);

public:
    static double pdf(double z);

    NormalMoment(double bounding);
};


inline double NormalMoment::pdf(double z) {
    const double x2 = z * z;
#if __cplusplus >= 202002L
    return std::exp(-0.5 * x2) / sqrt(2 * std::numbers::pi);
#else
    static const double pi = 3.14159265358979323846;
    return std::exp(-0.5 * x2) / sqrt(2 * pi);
#endif
}



inline NormalMoment::NormalMoment(double bounding)
#if __cplusplus >= 201103L
        : Moment(NormalMoment::init(bounding), bounding, 1 - std::erf(bounding/std::sqrt(2))) {
#else
        : Moment(NormalMoment::init(bounding), bounding, 1 - ::erf(bounding/std::sqrt(2.0))) {
#endif
}

inline std::vector<double> NormalMoment::init(double bounding) {
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

    std::vector<double> sMoment(n << 1, 0);
    for (size_t i = 0; i < n; ++i)
        sMoment[i << 1] = sTerm[i];
    for (size_t j = 2; j < sMoment.size(); ++j) {
        for (size_t i = 0; i < n; ++i) {
            sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2;
            const double prev = sMoment[i << 1];
            sMoment[i << 1] += sTerm[i];
            if (prev == sMoment[i << 1]) {
                n = i;
                break;
            }
        }
        if (n <= 0)
            break;
    }
    // Normalize per Formula (2.2): divide by ∫ρ dz = erf(κ/√2) so ζ(0,κ)=1.
#if __cplusplus >= 201103L
    const double normFactor = std::erf(bounding / std::sqrt(2));
#else
    const double normFactor = ::erf(bounding / std::sqrt(2.0));
#endif
    for (size_t i = 0; i < sMoment.size(); ++i) {
        sMoment[i] /= normFactor;
    }
    return sMoment;
}



} // namespace var_dbl
#endif // __Moment_h__
