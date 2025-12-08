/*
The ulp function
*/
#include <cmath>
#include <concepts>
#include <limits>

#ifndef __ulp_h__
#define __ulp_h__
namespace var_dbl 
{

/*
 * @return  the rounding error to hold x as a double
 */
template <typename T> requires std::integral<T> 
double ulp(T x) 
{
    static_assert(sizeof(unsigned long long) * 8 > std::numeric_limits<double>::digits);
    constexpr const unsigned long long max_significand 
        = (1LL << std::numeric_limits<double>::digits) - 1;
    unsigned long long val = (x >= 0)? x : -x;
    double round_error = 0;
    bool posi = true;
    while (val > max_significand) {
        round_error *= 0.5;
        if ((val & 1) == 1) {
            if (posi) {
                round_error += 1;
                posi = false;
            } else {
                round_error -= 1;
                posi = true;
            }
        }
        val >>= 1;
    }
    return round_error;   
}

/*
 * @return  units in last place of T type.
 */
template <typename T> requires std::floating_point<T> 
inline double ulp(T x) 
{
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<T>::infinity()) - x;
    else 
        return x - std::nexttoward(x, -std::numeric_limits<T>::infinity());   
}

/*
 * @return  If all least significant precise_tail_bits is 0, return 0;
            Otherwise, return units in last place of T type.
 */
template <typename T> requires std::floating_point<T> 
double ulp(T value, unsigned precise_tail_bits)
{
    if (precise_tail_bits > 0) {
        unsigned long mask = (1 << precise_tail_bits) - 1;
        if (sizeof(T) == sizeof(unsigned long long)) {
            unsigned long long u;
            memcpy( &u, &value, sizeof(u));
            if ((u & mask) == 0)
                return 0;
        } else if (sizeof(T) == sizeof(unsigned long)) {
            unsigned long u;
            memcpy( &u, &value, sizeof(u));
            if ((u & mask) == 0)
                return 0;
        } else if (sizeof(T) == sizeof(unsigned int)) {
            unsigned int u;
            memcpy( &u, &value, sizeof(u));
            if ((u & mask) == 0)
                return 0;
        } else if (sizeof(T) == sizeof(unsigned short)) {
            unsigned short u;
            memcpy( &u, &value, sizeof(u));
            if ((u & mask) == 0)
                return 0;
        }
    }
    return ulp(value);
}


/*
double ulp(double x)
{
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<double>::infinity()) - x;
    else 
        return x - std::nexttoward(x, -std::numeric_limits<double>::infinity());
}

double ulp(float x)
{
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<float>::infinity()) - x;
    else 
        return x - std::nexttoward(x, -std::numeric_limits<float>::infinity());
}
*/


} // namespace var_dbl
#endif // __ulp_h__