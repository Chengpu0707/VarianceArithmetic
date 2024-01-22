#include <cmath>
#include <exception>
#include <limits>
#include <format>   // not supported in gcc 13.2.0
#include <sstream>
#include <string>

#ifndef __ValDbl_h__
#define __ValDbl_h__

namespace var_dbl 
{

class ValueError : public std::exception 
{
public:
    const double value;
    const std::string what;

    explicit ValueError(double value, const std::string what) : 
        value(value), what(what)
    {
    }
};


class UncertaintyError : public std::exception 
{
public:
    const double value;
    const double variance;
    const std::string what;

    explicit UncertaintyError(double value, double variance, const std::string what) : 
        value(value), variance(variance), what(what)
    {
    }
};

struct VarDbl {     // a class which is cheap to copy
public:
    // assume uniform distribution within ulp()
    constexpr static const double DEVIATION_OF_LSB = 1.0 / sqrt(3);
    static double ulp(double d);


private:    
    double _value = 0;
    double _variance = 0;

    void init(double value, double variance, const std::string what) 
    {
        if (!std::isfinite(value))
            throw ValueError(value, what);
        if (!std::isfinite(variance))
            throw UncertaintyError(value, variance, what);
        _value = value;
        _variance = variance;
    }


public:
    double value() const { return _value; }
    double variance() const { return _variance; }
    double uncertainty() const { return sqrt(variance()); }

    // constructors
    VarDbl();
    VarDbl(const VarDbl& other);
    VarDbl(double value, double uncertainty);
        // uncertainty is limited between sqrt(std::numeric_limits<double>::mim()) and sqrt(std::numeric_limits<double>::max())
    // conversion constructors
    VarDbl(double value);
        // assume ulp as uncertainty
    VarDbl(float value);
        // assume ulp as uncertainty
    VarDbl(long long value) noexcept;
        // may loss resolution
    VarDbl(long value) noexcept : VarDbl((long long) value) {}
    VarDbl(int value) noexcept : VarDbl((long long) value) {}

    // i/o
    std::string to_string() const;
    friend std::ostream & operator <<(std::ostream& out, const VarDbl& v);
    friend std::istream & operator >>(std::istream& in, VarDbl& v);

    // +/-
    void negate() { _value = - _value; }
    VarDbl operator-() const;
    VarDbl operator+(VarDbl other) const;
    VarDbl operator-(VarDbl other) const;
    VarDbl operator+=(VarDbl other);
    VarDbl operator-=(VarDbl other);
    template<typename T> friend VarDbl operator+(T first, VarDbl second);
    template<typename T> friend VarDbl operator-(T first, VarDbl second);

    // *
    VarDbl operator*(VarDbl other) const;
    VarDbl operator*=(VarDbl other);
    template<typename T> friend VarDbl operator*(T first, VarDbl second);

};


inline double VarDbl::ulp(double x)
{
    if (x > 0)
        return (std::nexttoward(x, std::numeric_limits<double>::infinity()) - x) * DEVIATION_OF_LSB;
    else 
        return (x - std::nexttoward(x, -std::numeric_limits<double>::infinity())) * DEVIATION_OF_LSB;
}

inline VarDbl::VarDbl() {
    _value = 0;
    _variance = 0;
}

inline VarDbl::VarDbl(const VarDbl& other) 
{
    _value = other._value;
    _variance = other._variance;
}

inline VarDbl::VarDbl(double value, double uncertainty) 
{
    std::ostringstream ss;
    ss << "VarDbl(double " << value << ", double " << uncertainty << ")";
    init(value, uncertainty*uncertainty, ss.str());
}

inline VarDbl::VarDbl(double value) 
{
    std::ostringstream ss;
    ss << "VarDbl(double " << value << ")";
    const double uncertainty = ulp(value); 
    init(value, uncertainty*uncertainty, ss.str());
}

inline VarDbl::VarDbl(float value)
{
    std::ostringstream ss;
    ss << "VarDbl(float " << value << ")";
    const double uncertainty = DEVIATION_OF_LSB * ((value > 0)
        ? std::nexttoward(value, std::numeric_limits<float>::infinity()) - value
        : value - std::nexttoward(value, -std::numeric_limits<float>::infinity()));
    init(value, uncertainty*uncertainty, ss.str());
}

inline VarDbl::VarDbl(long long value) noexcept 
{
    if (value == (long long) ((double) value)) {
        _value = value;
        _variance = 0;
        return;
    }
    std::ostringstream ss;
    ss << "VarDbl(long " << value << ")";
    value = (double) value;
    const double uncertainty = ulp(value);
    init(value, uncertainty*uncertainty, ss.str());
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
    v._variance = uncertainty * uncertainty;
    return in;
}

inline VarDbl VarDbl::operator-() const 
{
    VarDbl ret(*this);
    ret.negate();
    return ret;
}

inline VarDbl VarDbl::operator+(VarDbl other) const
{
    other.init(_value + other._value, variance() + other.variance(), "+");
    return other;
}

inline VarDbl VarDbl::operator-(VarDbl other) const
{
    other.init(_value - other._value, variance() + other.variance(), "-");
    return other;
}

inline VarDbl VarDbl::operator+=(VarDbl other) 
{
    init(_value + other._value, variance() + other.variance(), "+=");
    return *this;
}

inline VarDbl VarDbl::operator-=(VarDbl other)
{
    init(_value - other._value, variance() + other.variance(), "-=");
    return *this;
}

template<typename T> 
inline VarDbl operator+(T first, VarDbl second)
{
    return second + first;
}

template<typename T> 
inline VarDbl operator-(T first, VarDbl second) 
{
    second -= first;
    second.negate();
    return second;
}

inline VarDbl VarDbl::operator*(VarDbl other) const
{
    other *= *this;
    return other;
}

inline VarDbl VarDbl::operator*=(VarDbl other) {
    const double variance = this->variance() * other.value() * other.value() +
                            other.variance() * value() * value() +
                            this->variance() * other.variance();
    init(value() * other.value(), variance, "*=");
    return *this;
}

template<typename T> 
inline VarDbl operator*(T first, VarDbl second) {
    return second * VarDbl(first);
}


} // namespace var_dbl
#endif // __ValDbl_h__
