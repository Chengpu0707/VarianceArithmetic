/*
The ulp function
*/

#include <cmath>
#include <limits>

#ifndef __ulp_h__
#define __ulp_h__
namespace var_dbl 
{

inline double ulp(double x)
{
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<double>::infinity()) - x;
    else 
        return x - std::nexttoward(x, -std::numeric_limits<double>::infinity());
}

} // namespace var_dbl
#endif // __ulp_h__