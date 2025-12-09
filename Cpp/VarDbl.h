/*
VarDbl is a class to implement variance arithmetic in C++

VarDbl.h is intended to be fully inline and deplored along without any making system.
It relies only on Momentum.h, which is also in the var_dbl namespace.
*/
#include <cmath>

#include "Test.h"
#include "ulp.h"

#ifndef __ValDbl_h__
#define __ValDbl_h__
namespace var_dbl 
{


class VarDbl { 
    constexpr static const double DEVIATION_OF_LSB = 1.0 / sqrt(3);
        // assume uniform distribution within ulp()
    constexpr static const double BINDING_FOR_EQUAL = 0.67448975;
        // z for 50% probability of equal
    constexpr static const long PRECISE_SIGNIFICAND_TAIL_BITS = 23;
        // If the last 23 bits of significand are zero, the value is 2's fractional.
    constexpr static const long long DOUBLE_MAX_SIGNIFICAND = (1LL << std::numeric_limits<double>::digits) - 1;
        // max significand for double

    double _value = 0;
    double _uncertainty = 0;

    void init(double value, double uncertainty, const std::string what);

public:
    template <typename T> requires std::floating_point<T> static double ulp(T value);

    double value() const { return _value; }
    double uncertainty() const { return _uncertainty; }
    double variance() const { return _uncertainty * _uncertainty; }

    // constructors
    VarDbl();
    VarDbl(const VarDbl& other);
    VarDbl(double value, double uncertainty);
    
    // conversion constructors
    template <typename T> requires std::floating_point<T> VarDbl(T value);
    template <typename T> requires std::integral<T> VarDbl(T value);

    // i/o
    std::string to_string() const;
    friend std::ostream & operator <<(std::ostream& out, const VarDbl& v);
    friend std::istream & operator >>(std::istream& in, VarDbl& v);

    // +
    VarDbl operator+(const VarDbl& other) const;
    VarDbl operator+=(const VarDbl& other);
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator+(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator+=(const T& other);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        VarDbl operator+(const T& first, const VarDbl& second);

    // -
    void negate() { _value = - _value; }
    VarDbl operator-() const;
    VarDbl operator-(const VarDbl& other) const;
    VarDbl operator-=(const VarDbl& other);
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator-(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator-=(const T& other);
    template<typename T>  requires std::floating_point<T> || std::integral<T> friend 
        VarDbl operator-(const T& first, const VarDbl& second);

    // *
    VarDbl operator*(const VarDbl& other) const;
    VarDbl operator*=(const VarDbl& other);
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator*(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator*=(const T& other);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        VarDbl operator*(const T& first, const VarDbl& second);

    // /
    VarDbl operator/(const VarDbl& other) const;
    VarDbl operator/=(const VarDbl& other);
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        VarDbl operator/(const T& other) const;
    template<typename T>  requires std::floating_point<T> || std::integral<T> 
        VarDbl operator/=(const T& other);
    template<typename T>  requires std::floating_point<T> || std::integral<T> 
        friend VarDbl operator/(const T& first, const VarDbl& second);

    // compare
    bool operator==(const VarDbl& other) const;
    bool operator!=(const VarDbl& other) const;
    bool operator<(const VarDbl& other) const;
    bool operator>(const VarDbl& other) const;
    bool operator<=(const VarDbl& other) const;
    bool operator>=(const VarDbl& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        bool operator==(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        bool operator!=(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        bool operator<(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> 
        bool operator>(const T& other) const;
    template<typename T>  requires std::floating_point<T> || std::integral<T> 
        bool operator<=(const T& other) const;
    template<typename T>  requires std::floating_point<T> || std::integral<T> 
        bool operator>=(const T& other) const;
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        bool operator==(const T& first, const VarDbl& second);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        bool operator!=(const T& first, const VarDbl& second);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        bool operator<(const T& first, const VarDbl& second);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        bool operator>(const T& first, const VarDbl& second);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        bool operator<=(const T& first, const VarDbl& second);
    template<typename T> requires std::floating_point<T> || std::integral<T> friend 
        bool operator>=(const T& first, const VarDbl& second);

    static void assertEqual(const VarDbl& var, double value, double uncertainty, const std::string& msg = "",
                     double valueDelta = 0, double uncertaintyDelta = 0);

};




struct InitException : public std::runtime_error
{
    const double value;
    const double uncertainty;

    explicit InitException(double value, double uncertainty, const std::string what) : 
        runtime_error(what), value(value), uncertainty(uncertainty)
    {}
};


inline void VarDbl::init(double value, double uncertainty, const std::string what) {
    if (!std::isfinite(value) || !std::isfinite(uncertainty))
        throw InitException(value, uncertainty, what);
    _value = value;
    _uncertainty = abs(uncertainty);
}


inline VarDbl::VarDbl() {
    _value = 0;
    _uncertainty = 0;
}

inline VarDbl::VarDbl(const VarDbl& other) 
{
    _value = other._value;
    _uncertainty = other._uncertainty;
}

inline VarDbl::VarDbl(double value, double uncertainty) 
{
    std::ostringstream ss;
    ss << "VarDbl( " << value << ", " << uncertainty << ")";
    init(value, uncertainty, ss.str());
}

template <typename T> requires std::floating_point<T> 
inline double VarDbl::ulp(T value) {
    return var_dbl::ulp(value, VarDbl::PRECISE_SIGNIFICAND_TAIL_BITS) * VarDbl::DEVIATION_OF_LSB;
}

template <typename T> requires std::floating_point<T> 
inline VarDbl::VarDbl(T value) {
    std::ostringstream ss;
    ss << "VarDbl(fp " << value << ")";
    const double uncertainty = VarDbl::ulp(value);
    init(value, uncertainty, ss.str());
}

template <typename T> requires std::integral<T> 
inline VarDbl::VarDbl(T value) {
    std::ostringstream ss;
    ss << "VarDbl(i " << value << ")";
    const double uncertainty = var_dbl::ulp(value);
    init(value, uncertainty, ss.str());
}

inline std::string VarDbl::to_string() const {
    std::ostringstream os;
    os << std::scientific << value() << '~' << uncertainty();
    return os.str();
}

inline std::ostream & operator <<(std::ostream& out, const VarDbl& v) 
{
    out << std::scientific << v.value() << '~' << v.uncertainty();
    return out;
}

inline std::istream & operator >>(std::istream& in, VarDbl& v) 
{
    double value, uncertainty;
    char sep;
    in >> value >> sep >> uncertainty;
    if (sep != '~') {
        std::ostringstream ss;
        ss << "Invalid separator for " << value << sep << uncertainty;
        throw std::invalid_argument(ss.str());
    }
    v._value = value;
    v._uncertainty = uncertainty;
    return in;
}


inline VarDbl VarDbl::operator-() const 
{
    VarDbl ret(*this);
    ret.negate();
    return ret;
}

inline VarDbl VarDbl::operator+=(const VarDbl& other) 
{
    const double uncertainty = (this->uncertainty() ==0)? other.uncertainty() :
            (other.uncertainty() == 0)? this->uncertainty() :
            std::sqrt(variance() + other.variance());
    init(_value + other._value, uncertainty, "+=");
    return *this;
}

inline VarDbl VarDbl::operator-=(const VarDbl& other)
{
    const double uncertainty = (this->uncertainty() ==0)? other.uncertainty() :
            (other.uncertainty() == 0)? this->uncertainty() :
            std::sqrt(variance() + other.variance());
    init(_value - other._value, uncertainty, "-=");
    return *this;
}

inline VarDbl VarDbl::operator*=(const VarDbl& other) 
{
    const double uncertainty = std::sqrt(
            this->variance() * other.value() * other.value() +
            other.variance() * value() * value() +
            this->variance() * other.variance());
    if (uncertainty <= 0) {
        const long long val = ((long long) value()) * ((long long) other.value());
        if (VarDbl::DOUBLE_MAX_SIGNIFICAND < std::abs(val) &&
            value() < VarDbl::DOUBLE_MAX_SIGNIFICAND && other.value() < VarDbl::DOUBLE_MAX_SIGNIFICAND) {
            const VarDbl v(val);
            this->_value = v._value;
            this->_uncertainty = v._uncertainty;
        } else
            init(value() * other.value(), uncertainty, "*=");
    } else
        init(value() * other.value(), uncertainty, "*=");
    return *this;
}

// VarDbl VarDbl::operator/=(const VarDbl& other) in Taylor.h


inline VarDbl VarDbl::operator+(const VarDbl& other) const
{
    VarDbl res(*this);
    res += other;
    return res;
}

inline VarDbl VarDbl::operator-(const VarDbl& other) const
{
    VarDbl res(*this);
    res -= other;
    return res;
}

inline VarDbl VarDbl::operator*(const VarDbl& other) const
{
    VarDbl res(*this);
    res *= other;
    return res;
}

inline VarDbl VarDbl::operator/(const VarDbl& other) const
{
    VarDbl res(*this);
    res /= other;
    return res;
}


template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator+(const T& other) const
{
    return *this + VarDbl(other);
}

template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator-(const T& other) const
{
    return *this - VarDbl(other);
}

template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator*(const T& other) const
{
    return *this * VarDbl(other);
}

template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator/(const T& other) const
{
    return *this / VarDbl(other);
}


template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator+=(const T& other)
{
    return *this += VarDbl(other);
}

template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator-=(const T& other)
{
    return *this -= VarDbl(other);
}

template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator*=(const T& other)
{
    return *this *= VarDbl(other);
}

template <typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl VarDbl::operator/=(const T& other)
{
    return *this /= VarDbl(other);
}


template<typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl operator+(const T& first, const VarDbl& second)
{
    return VarDbl(first) + second;
}

template<typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl operator-(const T& first, const VarDbl& second) 
{
    return VarDbl(first) - second;
}

template<typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl operator*(const T& first, const VarDbl& second) {
    return VarDbl(first) * second;
}

template<typename T> requires std::floating_point<T> || std::integral<T>  
inline VarDbl operator/(const T& first, const VarDbl& second) 
{
    return VarDbl(first) / second;
}


inline bool VarDbl::operator==(const VarDbl& other) const
{
    VarDbl res(*this);    
    res-= other;
    return std::abs(res.value()) <= (VarDbl::BINDING_FOR_EQUAL * res.uncertainty());
}

inline bool VarDbl::operator!=(const VarDbl& other) const
{
    return !(*this == other);
}

inline bool VarDbl::operator<(const VarDbl& other) const
{
    if (*this == other)
        return false;
    return this->value() < other.value();
}

inline bool VarDbl::operator>(const VarDbl& other) const
{
    if (*this == other)
        return false;
    return this->value() > other.value();
}

inline bool VarDbl::operator<=(const VarDbl& other) const
{
    if (*this == other)
        return true;
    return this->value() < other.value();
}

inline bool VarDbl::operator>=(const VarDbl& other) const
{
    if (*this == other)
        return true;
    return this->value() > other.value();
}


template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool VarDbl::operator==(const T& other) const
{
    return *this == VarDbl(other);
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool VarDbl::operator!=(const T& other) const
{
    return *this != VarDbl(other);
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool VarDbl::operator<(const T& other) const
{
    return *this < VarDbl(other);
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool VarDbl::operator>(const T& other) const
{
    return *this > VarDbl(other);
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool VarDbl::operator>=(const T& other) const
{
    return *this >= VarDbl(other);
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool VarDbl::operator<=(const T& other) const
{
    return *this <= VarDbl(other);
}


template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool operator==(const T& first, const VarDbl& second)
{
    return VarDbl(first) == second;
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool operator!=(const T& first, const VarDbl& second)
{
    return VarDbl(first) != second;
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool operator<(const T& first, const VarDbl& second)
{
    return VarDbl(first) < second;
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool operator>(const T& first, const VarDbl& second)
{
    return VarDbl(first) > second;
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool operator<=(const T& first, const VarDbl& second)
{
    return VarDbl(first) <= second;
}

template<typename T> requires std::floating_point<T> || std::integral<T> 
inline bool operator>=(const T& first, const VarDbl& second)
{
    return VarDbl(first) >= second; 
}


inline void VarDbl::assertEqual(const VarDbl& var, double value, double uncertainty, 
            const std::string& msg,
            double deltaValue, double deltaUncertainty) {
    std::ostringstream oss;
    if (!msg.empty())
        oss << msg << ": ";
    oss << var.value() << "~" << var.uncertainty() << " vs " << value << "~" << uncertainty;
    if (deltaValue == 0)
        deltaValue = std::max(ulp(var.value()), ulp(value));
    test::assertAlmostEqual(var.value(), value, deltaValue, oss.str());
    if (deltaUncertainty == 0)
        deltaUncertainty = std::max(ulp(var.uncertainty()), ulp(uncertainty));
    test::assertAlmostEqual(var.uncertainty(), uncertainty, deltaUncertainty, oss.str());
}



} // namespace var_dbl
#endif // __ValDbl_h__
