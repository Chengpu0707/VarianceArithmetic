#include <cmath>
#include <exception>
#include <format>   // not supported in gcc 13.2.0
#include <sstream>
#include <string>

#include "Test.h"

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


class VarianceError : public std::exception 
{
public:
    const double value;
    const double variance;
    const std::string what;

    explicit VarianceError(double value, double variance, const std::string what) : 
        value(value), variance(variance), what(what)
    {
    }
};

struct VarDbl {     // a plain-old-type which is easy to copy
private:    
    double _value = 0;
    double _variance = 0;

    void init(double value, double variance, const std::string what) 
    {
        if (!std::isfinite(value))
            throw ValueError(value, what);
        if (!std::isfinite(variance))
            throw VarianceError(value, variance, what);
        _value = value;
        _variance = variance;
    }

    // assume uniform distribution within ulp()
    constexpr static const double DEVIATION_OF_LSB = 1.0 / sqrt(3);

public:
    double value() const { return _value; }
    double uncertainty() const { return sqrt(_variance); }

    // constructors
    VarDbl();
    VarDbl(const VarDbl& other);
    VarDbl(double value, double uncertainty);
        // uncertainty is limited between sqrt(std::numeric_limits<double>::mim()) and sqrt(std::numeric_limits<double>::max())
    VarDbl(double value);
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
    VarDbl operator+(const VarDbl other) const;
    VarDbl operator-(const VarDbl other) const;
    VarDbl operator+=(const VarDbl other);
    VarDbl operator-=(const VarDbl other);


};

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
    const double uncertainty = Test::ulp(value)*DEVIATION_OF_LSB; 
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
    const double uncertainty = Test::ulp(value)*DEVIATION_OF_LSB;
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

inline VarDbl VarDbl::operator+(const VarDbl other) const
{
    VarDbl ret(*this);
    ret += other;
    return ret;
}

inline VarDbl VarDbl::operator-(const VarDbl other) const
{
    VarDbl ret(*this);
    ret -= other;
    return ret;
}

inline VarDbl VarDbl::operator+=(const VarDbl other) 
{
    init(_value + other._value, _variance + other._variance, "+=");
    return *this;
}

inline VarDbl VarDbl::operator-=(const VarDbl other)
{
    init(_value - other._value, _variance + other._variance, "-=");
    return *this;
}



} // namespace var_dbl
#endif // __ValDbl_h__
