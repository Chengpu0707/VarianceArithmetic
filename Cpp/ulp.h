/*
The ulp function
*/
#if __cplusplus >= 202002L
#include <bit>
#include <concepts>
#endif
#include <cmath>
#include <cstring>
#include <limits>
#if __cplusplus >= 201103L
#include <type_traits>
#endif

#ifndef __ulp_h__
#define __ulp_h__
namespace var_dbl
{

#if __cplusplus < 201103L
template<bool B, typename T> struct _enable_if_c03 {};
template<typename T> struct _enable_if_c03<true, T> { typedef T type; };
template<typename T> struct _is_floating_point_c03 { static const bool value = false; };
template<> struct _is_floating_point_c03<float> { static const bool value = true; };
template<> struct _is_floating_point_c03<double> { static const bool value = true; };
template<> struct _is_floating_point_c03<long double> { static const bool value = true; };
template<typename T> struct _is_integral_c03 { static const bool value = false; };
template<> struct _is_integral_c03<bool> { static const bool value = true; };
template<> struct _is_integral_c03<char> { static const bool value = true; };
template<> struct _is_integral_c03<signed char> { static const bool value = true; };
template<> struct _is_integral_c03<unsigned char> { static const bool value = true; };
template<> struct _is_integral_c03<short> { static const bool value = true; };
template<> struct _is_integral_c03<unsigned short> { static const bool value = true; };
template<> struct _is_integral_c03<int> { static const bool value = true; };
template<> struct _is_integral_c03<unsigned int> { static const bool value = true; };
template<> struct _is_integral_c03<long> { static const bool value = true; };
template<> struct _is_integral_c03<unsigned long> { static const bool value = true; };
template<> struct _is_integral_c03<long long> { static const bool value = true; };
template<> struct _is_integral_c03<unsigned long long> { static const bool value = true; };
template<> struct _is_integral_c03<wchar_t> { static const bool value = true; };
struct _C03FloatTag {};
struct _C03IntTag {};
#endif

/*
 * @return  the rounding error to hold x as a double
 */
#if __cplusplus >= 202002L
template <typename T> requires std::integral<T>
double ulp(T x)
#elif __cplusplus >= 201103L
template <typename T>
typename std::enable_if<std::is_integral<T>::value, double>::type
ulp(T x)
#else
template <typename T>
typename _enable_if_c03<_is_integral_c03<T>::value, double>::type
ulp(T x)
#endif
{
#if __cplusplus >= 201103L
    static_assert(sizeof(unsigned long long) * 8 > std::numeric_limits<double>::digits);
    constexpr const unsigned long long max_significand
        = (1LL << std::numeric_limits<double>::digits) - 1;
#else
    const unsigned long long max_significand
        = (1ULL << std::numeric_limits<double>::digits) - 1;
#endif
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
        // round_error *= 0.5;
    }
    return round_error;
}

/*
 * @return  units in last place of T type.
 */
#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
inline double ulp(T x)
#elif __cplusplus >= 201103L
template <typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, double>::type
ulp(T x)
#else
template <typename T>
inline typename _enable_if_c03<_is_floating_point_c03<T>::value, double>::type
ulp(T x)
#endif
{
#if __cplusplus >= 201103L
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<T>::infinity()) - x;
    else
        return x - std::nexttoward(x, -std::numeric_limits<T>::infinity());
#else
    if (x > 0)
        return ::nexttoward(x, std::numeric_limits<T>::infinity()) - x;
    else
        return x - ::nexttoward(x, -std::numeric_limits<T>::infinity());
#endif
}

/*
 * @return  If all least significant precise_tail_bits is 0, return 0;
            Otherwise, return units in last place of T type.
 */
#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
double ulp(T value, unsigned precise_tail_bits)
#elif __cplusplus >= 201103L
template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, double>::type
ulp(T value, unsigned precise_tail_bits)
#else
template <typename T>
typename _enable_if_c03<_is_floating_point_c03<T>::value, double>::type
ulp(T value, unsigned precise_tail_bits)
#endif
{
    if (precise_tail_bits > 0) {
        unsigned long mask = (1 << precise_tail_bits) - 1;
#if __cplusplus >= 202002L
        if constexpr (sizeof(T) == sizeof(unsigned long long)) {
            if ((std::bit_cast<unsigned long long>(value) & mask) == 0)
                return 0;
        } else if constexpr (sizeof(T) == sizeof(unsigned long)) {
            if ((std::bit_cast<unsigned long>(value) & mask) == 0)
                return 0;
        } else if constexpr (sizeof(T) == sizeof(unsigned int)) {
            if ((std::bit_cast<unsigned int>(value) & mask) == 0)
                return 0;
        } else if constexpr (sizeof(T) == sizeof(unsigned short)) {
            if ((std::bit_cast<unsigned short>(value) & mask) == 0)
                return 0;
        }
#else
        if (sizeof(T) == sizeof(unsigned long long)) {
            unsigned long long _bits; std::memcpy(&_bits, &value, sizeof(_bits));
            if ((_bits & mask) == 0) return 0;
        } else if (sizeof(T) == sizeof(unsigned long)) {
            unsigned long _bits; std::memcpy(&_bits, &value, sizeof(_bits));
            if ((_bits & mask) == 0) return 0;
        } else if (sizeof(T) == sizeof(unsigned int)) {
            unsigned int _bits; std::memcpy(&_bits, &value, sizeof(_bits));
            if ((_bits & mask) == 0) return 0;
        } else if (sizeof(T) == sizeof(unsigned short)) {
            unsigned short _bits; std::memcpy(&_bits, &value, sizeof(_bits));
            if ((_bits & mask) == 0) return 0;
        }
#endif
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
