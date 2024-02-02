/*
Momentum is needed for Taylor expansion in VarDbl.  
It is a constexpr class.

Unlike in the python version, it is in double only, because the uncertainty is close to LSB of significand.
Also, the double format is more friendly to be incorporated into the C++ headers.
*/

#include <cmath>

#ifndef __Momentum_h__
#define __Momentum_h__
namespace var_dbl 
{

template<int maxOrder, int binding>
class Momentum {
public:
    double factor(int n) const 
    {
        if ((n % 2) == 1)
            return 0;
        n /= 2;
        if (n > maxOrder)
            return 0;
        return _sFactor[n];
    }

    constexpr Momentum() : _sFactor() {
        constexpr const int divid = 32;
        constexpr const double divid2 = divid * divid;
        constexpr const int limit = divid * binding;
        for (int j = 0; j < maxOrder; ++j)
            _sFactor[j] = 0;
        for (int i = -limit ; i <= limit; ++i) {
            const double x2 = i*i / divid2;
            const double pdf = 1.0/sqrt(2*M_PI) * exp(- x2 * 0.5) / divid;
            double sq = 1;
            for (int j = 0; j < maxOrder; ++j)
            {
                _sFactor[j] += pdf * sq;
                sq *= x2;
            }
        }
    }

private:
    double _sFactor[maxOrder];
};


} // namespace var_dbl
#endif // __Momentum_h__