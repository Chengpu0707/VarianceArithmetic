/*
VarDbl is a class to implement variance arithmetic in C++
It is cheap to copy or to create, to encourage type conversion by copying. 

VarDbl.h is intended to be fully inline and deplored along without any making system.
It relies only on Momentum.h, which is also in the var_dbl namespace.
*/
#include <cmath>
#include <exception>
#include <format>   // not supported in gcc 13.2.0
#include <sstream>
#include <string>

#include "Momentum.h"
#include "ulp.h"

#ifndef __ValDbl_h__
#define __ValDbl_h__
namespace var_dbl 
{

struct ValueException : public std::runtime_error
{
    const double value;

    explicit ValueException(double value, const std::string what) : 
        runtime_error(what), value(value)
    {}
};


struct UncertaintyException : public std::runtime_error 
{
    const double value;
    const double variance;

    explicit UncertaintyException(double value, double variance, const std::string what) : 
        runtime_error(what), value(value), variance(variance)
    {}
};

struct VarDbl { 
    constexpr static const double DEVIATION_OF_LSB = 1.0 / sqrt(3);
        // assume uniform distribution within ulp()
    constexpr static const double BINDING_FOR_EQUAL = 0.67448975;
        // z for 50% probability of equal

    // Taylor expansion
    constexpr static const int MAX_ORDER_FOR_TAYLOR = 126;
    constexpr static const int BINDING_FOR_TAYLOR = 5;
    constexpr static const double VARIANCE_THRESHOLD = 1.0 /BINDING_FOR_TAYLOR /BINDING_FOR_TAYLOR;

private:    
    constexpr static const auto _momentum = Momentum<MAX_ORDER_FOR_TAYLOR, BINDING_FOR_TAYLOR>();
    double _value = 0;
    double _variance = 0;

    void init(double value, double variance, const std::string what) 
    {
        if (!std::isfinite(value))
            throw ValueException(value, what);
        if (!std::isfinite(variance))
            throw UncertaintyException(value, variance, what);
        _value = value;
        _variance = variance;
    }

    VarDbl taylor1d(std::string name, VarDbl sTaylor[], bool inPrec, bool outPrec) const;

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
        // may result in uncertainty
    VarDbl(long value) noexcept : VarDbl((long long) value) {}
    // always without uncertainty
    VarDbl(int value) noexcept : VarDbl((long long) value) {}
    VarDbl(short value) noexcept : VarDbl((long long) value) {}
    VarDbl(char value) noexcept : VarDbl((long long) value) {}

    // i/o
    std::string to_string() const;
    friend std::ostream & operator <<(std::ostream& out, const VarDbl& v);
    friend std::istream & operator >>(std::istream& in, VarDbl& v);

    // +
    VarDbl operator+(VarDbl other) const;
    VarDbl operator+=(VarDbl other);
    template<typename T> friend VarDbl operator+(T first, VarDbl second);

    // -
    void negate() { _value = - _value; }
    VarDbl operator-() const;
    VarDbl operator-(VarDbl other) const;
    VarDbl operator-=(VarDbl other);
    template<typename T> friend VarDbl operator-(T first, VarDbl second);

    // *
    VarDbl operator*(VarDbl other) const;
    VarDbl operator*=(VarDbl other);
    template<typename T> friend VarDbl operator*(T first, VarDbl second);

    // compare
    bool operator==(VarDbl other) const;
    bool operator!=(VarDbl other) const;
    bool operator<(VarDbl other) const;
    bool operator>(VarDbl other) const;
    bool operator<=(VarDbl other) const;
    bool operator>=(VarDbl other) const;
    template<typename T> friend bool operator==(T first, VarDbl second);
    template<typename T> friend bool operator!=(T first, VarDbl second);
    template<typename T> friend bool operator<(T first, VarDbl second);
    template<typename T> friend bool operator>(T first, VarDbl second);
    template<typename T> friend bool operator<=(T first, VarDbl second);
    template<typename T> friend bool operator>=(T first, VarDbl second);

    // Taylor expansion, may throw LossUncertaintyException
    VarDbl exp() const;
    VarDbl log() const;
    VarDbl sin() const;
    VarDbl pow(double exp) const;

};


struct LossUncertaintyException : public std::runtime_error 
{
    const VarDbl input;
    const bool inPrec, outPrec;
    const VarDbl value;
    const VarDbl variance;
    const size_t n;
    const VarDbl newValue;
    const VarDbl newVariance;

    explicit LossUncertaintyException(const std::string what, 
                const VarDbl& input, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance, 
                size_t n, const VarDbl& newValue, const VarDbl& newVariance) :
            runtime_error(what + " LossUncertaintyExceptio"), input(input), inPrec(inPrec), outPrec(outPrec),
            value(value), variance(variance), n(n), newValue(newValue), newVariance(newVariance)
    {}
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
    const double uncertainty = ulp(value) * DEVIATION_OF_LSB; 
    ss << "VarDbl(double " << value << "): uncertainty=" << uncertainty;
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
    const double uncertainty = ulp(value) * DEVIATION_OF_LSB;
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

inline VarDbl VarDbl::operator*=(VarDbl other) 
{
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



inline bool VarDbl::operator==(VarDbl other) const
{
    other-= *this;
    if (other.value() == 0)
        return true;
    const double uncertainty = other.uncertainty();
    if (uncertainty == 0)
        return false;
    return abs(other.value() / uncertainty) <= VarDbl::BINDING_FOR_EQUAL;
}

inline bool VarDbl::operator!=(VarDbl other) const
{
    return !(*this == other);
}

inline bool VarDbl::operator<(VarDbl other) const
{
    if (*this == other)
        return false;
    return this->value() < other.value();
}

inline bool VarDbl::operator>(VarDbl other) const
{
    if (*this == other)
        return false;
    return this->value() > other.value();
}

inline bool VarDbl::operator<=(VarDbl other) const
{
    if (*this == other)
        return true;
    return this->value() < other.value();
}

inline bool VarDbl::operator>=(VarDbl other) const
{
    if (*this == other)
        return true;
    return this->value() > other.value();
}

template<typename T>
inline bool operator==(T first, VarDbl second)
{
    return second == first;
}

template<typename T>
inline bool operator!=(T first, VarDbl second)
{
    return second != first;
}

template<typename T>
inline bool operator<(T first, VarDbl second)
{
    return second > first;
}

template<typename T>
inline bool operator>(T first, VarDbl second)
{
    return second < first;
}

template<typename T>
inline bool operator<=(T first, VarDbl second)
{
    return second >= first;
}

template<typename T>
inline bool operator>=(T first, VarDbl second)
{
    return second <= first; 
}


inline VarDbl VarDbl::taylor1d(std::string name, VarDbl s1dTaylor[], bool inPrec, bool outPrec) const
{
    VarDbl value = outPrec? VarDbl(1, 0) : s1dTaylor[0];
    VarDbl variance = VarDbl();
    VarDbl var = VarDbl(this->variance());
    if (inPrec)
        var *= VarDbl( 1.0 /this->value() /this->value() );
    VarDbl varn = VarDbl(var);
    for (size_t n = 2; n < MAX_ORDER_FOR_TAYLOR; n += 2) {
        VarDbl newValue = s1dTaylor[n] * _momentum.factor(n) * varn;
        VarDbl newVariance = VarDbl();
        for (size_t j = 1; j < n; ++j) {
            newVariance += s1dTaylor[j] * s1dTaylor[n - j] * _momentum.factor(n) * varn;
        }
        for (size_t j = 2; j < n; j += 2) {
            newVariance -= s1dTaylor[j] * _momentum.factor(j) * \
                        s1dTaylor[n - j] * _momentum.factor(n - j) * varn;
        }
        value += newValue;
        variance += newVariance;
        if (variance.variance() > variance.value() * VARIANCE_THRESHOLD) {
            throw LossUncertaintyException(name, *this, inPrec, outPrec,
                    value, variance, n, newValue, newVariance);
        }
        varn *= var;
    }
    if (outPrec) {
        value *= s1dTaylor[0];
        variance *= s1dTaylor[0] * s1dTaylor[0];
    }
    VarDbl res;
    res.init(value.value(), variance.value() + value.variance(), name);
    return res;
}

inline VarDbl VarDbl::exp() const
{
    VarDbl sTaylor[MAX_ORDER_FOR_TAYLOR];
    sTaylor[0] = VarDbl(std::exp(value()));
    VarDbl n = VarDbl(1, 0);
    for (size_t i = 1; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        n *= 1.0 / i;
        sTaylor[i] = n;
    }
    std::ostringstream os;
    os << *this << ".exp()";
    return taylor1d(os.str(), sTaylor, false, true);
}

inline VarDbl VarDbl::log() const 
{
    VarDbl sTaylor[MAX_ORDER_FOR_TAYLOR];
    sTaylor[0] = VarDbl(std::log(value()));
    for (size_t i = 1; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        sTaylor[i] = VarDbl((((i%2) == 1)? 1.0 : -1.0) / i);
    }
    std::ostringstream os;
    os << *this << ".log()";
    return taylor1d(os.str(), sTaylor, true, false);
}

inline VarDbl VarDbl::sin() const 
{
    VarDbl sTaylor[MAX_ORDER_FOR_TAYLOR];
    sTaylor[0] = VarDbl(std::sin(value()));
    VarDbl n = VarDbl(1, 0);
    for (size_t i = 1; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        n *= 1.0 / i;
        switch (i % 4) {
            case 0:
                sTaylor[i] = VarDbl(std::sin(value()) * n);
                break;
            case 1:
                sTaylor[i] = VarDbl(std::cos(value()) * n);  
                break;
            case 2:
                sTaylor[i] = VarDbl(-std::sin(value()) * n);
                break;
            case 3:
                sTaylor[i] = VarDbl(-std::cos(value()) * n);
                break;
        }
    }
    std::ostringstream os;
    os << *this << ".sin()";
    return taylor1d(os.str(), sTaylor, false, false);
}

inline VarDbl VarDbl::pow(double exp) const 
{
    
}

} // namespace var_dbl
#endif // __ValDbl_h__
