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
#include <vector>

#include "Momentum.h"
#include "ulp.h"

#ifndef __ValDbl_h__
#define __ValDbl_h__
namespace var_dbl 
{

struct InitException : public std::runtime_error
{
    const double value;
    const double variance;

    explicit InitException(double value, double variance, const std::string what) : 
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
    constexpr static const int TAYLOR_CHECK_ORDER = 20;
    constexpr static const double TAU = 7.18e-7;  
    constexpr static const long PRECISE_SIGNIFICAND_TAIL_BITS = 23;

    template <typename T> requires std::floating_point<T> static double ulp(T value);

private:    
    constexpr static const auto _momentum = Momentum<MAX_ORDER_FOR_TAYLOR, BINDING_FOR_TAYLOR>();
    double _value = 0;
    double _variance = 0;

    void init(double value, double variance, const std::string what) {
        if (!std::isfinite(value) || !std::isfinite(variance))
            throw InitException(value, variance, what);
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
    VarDbl(double value, double uncertainty, bool uncertaintyAsVariance = false);
    
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
    template<typename T> friend VarDbl operator+(T first, const VarDbl& second);

    // -
    void negate() { _value = - _value; }
    VarDbl operator-() const;
    VarDbl operator-(const VarDbl& other) const;
    VarDbl operator-=(const VarDbl& other);
    template<typename T> friend VarDbl operator-(T first, const VarDbl& second);

    // *
    VarDbl operator*(const VarDbl& other) const;
    VarDbl operator*=(const VarDbl& other);
    template<typename T> friend VarDbl operator*(T first, const VarDbl& second);

    // *
    VarDbl operator/(const VarDbl& other) const;
    VarDbl operator/=(const VarDbl& other);
    template<typename T> friend VarDbl operator/(T first, const VarDbl& second);

    // compare
    bool operator==(const VarDbl& other) const;
    bool operator!=(const VarDbl& other) const;
    bool operator<(const VarDbl& other) const;
    bool operator>(const VarDbl& other) const;
    bool operator<=(const VarDbl& other) const;
    bool operator>=(const VarDbl& other) const;
    template<typename T> friend bool operator==(T first, const VarDbl& second);
    template<typename T> friend bool operator!=(T first, const VarDbl& second);
    template<typename T> friend bool operator<(T first, const VarDbl& second);
    template<typename T> friend bool operator>(T first, const VarDbl& second);
    template<typename T> friend bool operator<=(T first, const VarDbl& second);
    template<typename T> friend bool operator>=(T first, const VarDbl& second);

    // Taylor expansion, may throw NotReliableException
    VarDbl exp() const;
    VarDbl log() const;
    VarDbl sin() const;
    VarDbl pow(double exp) const;


    /*
     * 1d Taylor expansion.
     * 
     * @see     The paper for Variance Arithmetic on Taylor expansion convergence.
     * 
     * @return  The result of taylor expansion with this as input.
     *          If the input is precise, the uncertainty is the floating-point uncertainty of the result.
     * 
     * @param name          The name of the Taylor expansion, for exception logging.
     * @param s1dTaylor     The Taylor expansion coefficent, with f(x) as s1dTaylor[0]. It should already contains /n!.
     * @param inPrec        If to expand by input precision
     * @param outPrec       If the variance result needs to be multiplied by s1dTaylor[0].
     * @param dumpPath                      If to dump the expansion to a file
     * @param enableStabilityTruncation     If to truncate when the expansion becomes stable.
     * 
     * @exception InitException         If any item in s1dTaylor is not finite.
     * @exception DivergentException    If the result is not finite.
     * @exception NotReliableException  If the uncertainty of the variance is too large for its value. 
     * @exceptopm NotMonotonicException If the result variance does not decrease monotonically. 
     * @exceptopm NotStableException    If after maximal order expansion, the expansion is still not stable.       
     */
    VarDbl taylor1d(std::string name, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                    bool enableStabilityTruncation = true) const;

    /*
     * 1d polynominal
     * 
     * @see     The paper for Variance Arithmetic on Taylor expansion convergence.
     * 
     * @return  The result of taylor expansion with this as input.
     * 
     * @param sCoeff                    Polynominal coefficients
     * @param dumpPath                  If to dump the expansion to a file
     * 
     * @exception InitException         If any item in sCoeff is not finite.
     * @exception DivergentException    If the result is not finite.
     * @exception NotReliableException  If the uncertainty of the variance is too large for its value. 
     * 
     * @raise IllegalArgumentException  If the length of sCoeff is too long, or if sCoeff contains null
     */
    VarDbl polynominal(const std::vector<VarDbl>& sCoeff) const;

};


struct NotReliableException : public std::runtime_error 
{
    const VarDbl input;
    const std::vector<VarDbl> s1dTaylor;
    const bool inPrec, outPrec;
    const VarDbl value;
    const VarDbl variance;
    const size_t n;
    const VarDbl newValue;
    const VarDbl newVariance;

    explicit NotReliableException(const std::string what, 
                const VarDbl& input, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance, 
                size_t n, const VarDbl& newValue, const VarDbl& newVariance) :
            runtime_error(what + " NotReliableException"), input(input), 
            s1dTaylor(s1dTaylor), inPrec(inPrec), outPrec(outPrec),
            value(value), variance(variance), n(n), newValue(newValue), newVariance(newVariance)
    {}
};


struct DivergentException : public std::runtime_error 
{
    const VarDbl input;
    const bool inPrec, outPrec;
    const std::vector<VarDbl> s1dTaylor;
    const VarDbl value;
    const VarDbl variance;
    const size_t n;
    const VarDbl newValue;
    const VarDbl newVariance;

    explicit DivergentException(const std::string what, 
                const VarDbl& input, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance, 
                size_t n, const VarDbl& newValue, const VarDbl& newVariance) :
            runtime_error(what + " DivergentException"), input(input), 
            s1dTaylor(s1dTaylor), inPrec(inPrec), outPrec(outPrec),
            value(value), variance(variance), n(n), newValue(newValue), newVariance(newVariance)
    {}
};


struct NotMonotonicException : public std::runtime_error 
{
    const VarDbl input;
    const bool inPrec, outPrec;
    const std::vector<VarDbl> s1dTaylor;
    const VarDbl value;
    const VarDbl variance;
    const size_t n;
    const VarDbl newValue;
    const VarDbl newVariance;

    explicit NotMonotonicException(const std::string what, 
                const VarDbl& input, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance, 
                size_t n, const VarDbl& newValue, const VarDbl& newVariance) :
            runtime_error(what + " NotMonotonicException"), input(input), 
            s1dTaylor(s1dTaylor), inPrec(inPrec), outPrec(outPrec),
            value(value), variance(variance), n(n), newValue(newValue), newVariance(newVariance)
    {}
};


struct NotStableException : public std::runtime_error 
{
    const VarDbl input;
    const bool inPrec, outPrec;
    const std::vector<VarDbl> s1dTaylor;
    const VarDbl value;
    const VarDbl variance;
    const size_t n;
    const VarDbl newValue;
    const VarDbl newVariance;

    explicit NotStableException(const std::string what, 
                const VarDbl& input, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance, 
                size_t n, const VarDbl& newValue, const VarDbl& newVariance) :
            runtime_error(what + " NotStableException"), input(input), 
            s1dTaylor(s1dTaylor), inPrec(inPrec), outPrec(outPrec),
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

inline VarDbl::VarDbl(double value, double uncertainty, bool uncertaintyAsVariance) 
{
    std::ostringstream ss;
    ss << "VarDbl( " << value << ", " << (uncertaintyAsVariance? "var " : "unc ")  << uncertainty << ")";
    init(value, uncertaintyAsVariance? uncertainty : uncertainty*uncertainty, ss.str());
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
    init(value, uncertainty*uncertainty, ss.str());
}

template <typename T> requires std::integral<T> 
inline VarDbl::VarDbl(T value) {
    std::ostringstream ss;
    ss << "VarDbl(i " << value << ")";
    const double uncertainty = var_dbl::ulp(value);
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

inline VarDbl VarDbl::operator+(const VarDbl& other) const
{
    VarDbl res;
    res.init(_value + other._value, variance() + other.variance(), "+");
    return res;
}

inline VarDbl VarDbl::operator-(const VarDbl& other) const
{
    VarDbl res;
    res.init(_value - other._value, variance() + other.variance(), "-");
    return res;
}

inline VarDbl VarDbl::operator+=(const VarDbl& other) 
{
    init(_value + other._value, variance() + other.variance(), "+=");
    return *this;
}

inline VarDbl VarDbl::operator-=(const VarDbl& other)
{
    init(_value - other._value, variance() + other.variance(), "-=");
    return *this;
}

template<typename T> 
inline VarDbl operator+(T first, const VarDbl& second)
{
    return second + first;
}

template<typename T> 
inline VarDbl operator-(T first, const VarDbl& second) 
{
    VarDbl res(first);
    res.init(res.value() - second.value(), res.variance() + second.variance(), "-");
    return res;
}

inline VarDbl VarDbl::operator*(const VarDbl& other) const
{
    VarDbl res(other);
    const double variance = this->variance() * other.value() * other.value() +
                            other.variance() * value() * value() +
                            this->variance() * other.variance();
    res.init(value() * other.value(), variance, "*");
    return res;
}

inline VarDbl VarDbl::operator*=(const VarDbl& other) 
{
    const double variance = this->variance() * other.value() * other.value() +
                            other.variance() * value() * value() +
                            this->variance() * other.variance();
    init(value() * other.value(), variance, "*=");
    return *this;
}

template<typename T> 
inline VarDbl operator*(T first, const VarDbl& second) {
    return second * VarDbl(first);
}



inline VarDbl VarDbl::operator/(const VarDbl& other) const
{
    VarDbl res(*this);    
    return res * other.pow(-1);
}

inline VarDbl VarDbl::operator/=(const VarDbl& other) 
{
    return *this *= other.pow(-1);
}

template<typename T> 
inline VarDbl operator/(T first, const VarDbl& second) {
    return VarDbl(first) * second.pow(-1);
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

template<typename T>
inline bool operator==(T first, const VarDbl& second)
{
    return second == first;
}

template<typename T>
inline bool operator!=(T first, const VarDbl& second)
{
    return second != first;
}

template<typename T>
inline bool operator<(T first, const VarDbl& second)
{
    return second > first;
}

template<typename T>
inline bool operator>(T first, const VarDbl& second)
{
    return second < first;
}

template<typename T>
inline bool operator<=(T first, const VarDbl& second)
{
    return second >= first;
}

template<typename T>
inline bool operator>=(T first, const VarDbl& second)
{
    return second <= first; 
}


inline VarDbl VarDbl::taylor1d(std::string name, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                               bool enableStabilityTruncation) const
{
    if (variance() == 0) {
        return VarDbl(s1dTaylor[0]);
    }
    VarDbl value = outPrec? VarDbl(1, 0) : s1dTaylor[0];
    VarDbl variance = VarDbl();
    VarDbl var = VarDbl(this->variance());
    if (inPrec)
        var *= VarDbl( 1.0 /this->value() /this->value() );
    VarDbl varn = VarDbl(var);
    VarDbl prevVariance;
    for (size_t n = 2; 
            (n < MAX_ORDER_FOR_TAYLOR) && (n < s1dTaylor.size()) && (varn.value() > 0); 
            n += 2, varn *= var) {
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
        const double unc = variance.value() * VarDbl::TAU * VarDbl::TAU;

        if ((! std::isfinite(value.value())) || (! std::isfinite(value.variance()))
                || (! std::isfinite(variance.value())) || (! std::isfinite(variance.variance()))) {
            throw DivergentException(name, *this, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance);
        }
        if (variance.variance() > variance.value() * VARIANCE_THRESHOLD) {
            throw NotReliableException(name, *this, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance);
        }
        if (n >= TAYLOR_CHECK_ORDER) {
            if (std::abs(prevVariance.value()) + unc < std::abs(newVariance.value())) {
                throw NotMonotonicException(name, *this, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance);
            }
            if (enableStabilityTruncation) {
                if ((std::abs(newVariance.value()) < unc) || (std::abs(newValue.value()) < var_dbl::ulp(value.value())))
                    break;
                if (n == (MAX_ORDER_FOR_TAYLOR - 1)) {
                    throw NotStableException(name, *this, s1dTaylor, inPrec, outPrec,
                            value, variance, n, newValue, newVariance);
                }
            }
        } 
        prevVariance = newVariance;
    }
    if (outPrec) {
        value *= s1dTaylor[0];
        variance *= s1dTaylor[0] * s1dTaylor[0];
    }
    VarDbl res;
    res.init(value.value(), variance.value() + value.variance(), name);
    return res;
}


inline VarDbl VarDbl::polynominal(const std::vector<VarDbl>& sCoeff) const
{
    if (sCoeff.size() >= VarDbl::MAX_ORDER_FOR_TAYLOR*2) {
        std::ostringstream oss;
        oss << "coefficient length " << sCoeff.size() << " is too large!";
        throw std::invalid_argument(oss.str());
    }
    const int exp = sCoeff.size() - 1;
    std::vector<VarDbl> s1dTaylor(2*exp + 1);
    s1dTaylor[0] = sCoeff[0];
    std::vector<double> sPow;
    sPow.push_back(1);
    sPow.push_back(value());
    long sTaylor[2*exp + 1];
    sTaylor[0] = 1;
    for (int j = 1; j < sCoeff.size(); ++j) {
        const VarDbl& coeff = sCoeff[j];
        if (!std::isfinite(coeff.value()) || !std::isfinite(coeff.variance())) {
            std::ostringstream oss;
            oss << "polynominal coefficent at " << j << "/" << sCoeff.size();
            throw InitException(coeff.value(), coeff.variance(), oss.str());
        }
        sTaylor[1] = j;
        for (int k = 2; k <= j; ++k) {
            sTaylor[k] = sTaylor[k - 1] * (j + 1 - k) / k;
            if (sPow.size() <= k)
                sPow.push_back(sPow.back() * value());
        }
        for (int k = 0; k <= j; ++k) {
            s1dTaylor[k] += coeff * sTaylor[k] * sPow[j - k];
        }
    }

    VarDbl value = s1dTaylor[0];
    VarDbl variance = VarDbl();
    const VarDbl var = VarDbl(this->variance());
    VarDbl varn = VarDbl(var);
    for (int n = 2; (n < MAX_ORDER_FOR_TAYLOR*2) && (0 < varn.value()) && (n < s1dTaylor.size()); 
            n += 2, varn *= var) {
        VarDbl newValue = (n >= s1dTaylor.size())? VarDbl() : varn * s1dTaylor[n] * _momentum.factor(n);
        VarDbl newVariance = VarDbl();
        for (int j = 1; j < n; ++j) {
            newVariance += s1dTaylor[j] * s1dTaylor[n - j] * _momentum.factor(n) * varn;
        }
        for (int j = 2; j < n; j += 2) {
            newVariance -= s1dTaylor[j] * _momentum.factor(j) * s1dTaylor[n - j] * _momentum.factor(n - j) * varn;
        }
        value += newValue;
        variance += newVariance;

        if ((! std::isfinite(value.value())) || (! std::isfinite(value.variance()))
                || (! std::isfinite(variance.value())) || (! std::isfinite(variance.variance()))) {
            throw DivergentException("polynominal DivergentException", *this, s1dTaylor, false, false,
                    value, variance, n, newValue, newVariance);
        }
        if (variance.variance() > variance.value() * VARIANCE_THRESHOLD) {
            throw NotReliableException("polynominal NotReliableException", *this, s1dTaylor, false, false,
                    value, variance, n, newValue, newVariance);
        }
    }
    if ((variance.value() + value.variance()) == 0)
        return VarDbl(value.value());
    VarDbl res;
    res.init(value.value(), variance.value() + value.variance(), "polynominal");
    return res;
}


inline VarDbl VarDbl::exp() const
{
    std::vector<VarDbl> sTaylor;
    sTaylor.reserve(MAX_ORDER_FOR_TAYLOR);
    sTaylor.push_back(VarDbl(std::exp(value())));
    VarDbl n(1);
    for (size_t i = 1; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        n *= 1.0 / i;
        sTaylor.push_back(VarDbl(n));
    }
    std::ostringstream os;
    os << "exp(" << *this << ")";
    return taylor1d(os.str(), sTaylor, false, true);
}

inline VarDbl VarDbl::log() const 
{
    std::vector<VarDbl> sTaylor;
    sTaylor.reserve(MAX_ORDER_FOR_TAYLOR);
    sTaylor.push_back(VarDbl(std::log(value())));
    for (size_t i = 1; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        sTaylor.push_back(VarDbl((((i%2) == 1)? 1.0 : -1.0) / i));
    }
    std::ostringstream os;
    os << "log(" << *this << ")";
    return taylor1d(os.str(), sTaylor, true, false);
}

inline VarDbl VarDbl::sin() const 
{
    std::vector<VarDbl> sTaylor;
    sTaylor.reserve(MAX_ORDER_FOR_TAYLOR);
    sTaylor.push_back(VarDbl(std::sin(value())));
    VarDbl n = VarDbl(1, 0);
    for (size_t i = 1; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        n *= 1.0 / i;
        switch (i % 4) {
            case 0:
                sTaylor.push_back(VarDbl(std::sin(value()) * n));
                break;
            case 1:
                sTaylor.push_back(VarDbl(std::cos(value()) * n));  
                break;
            case 2:
                sTaylor.push_back(VarDbl(-std::sin(value()) * n));
                break;
            case 3:
                sTaylor.push_back(VarDbl(-std::cos(value()) * n));
                break;
        }
    }
    std::ostringstream os;
    os << "sin(" << *this << ")";
    return taylor1d(os.str(), sTaylor, false, false);
}

inline VarDbl VarDbl::pow(double exp) const 
{
    switch ((size_t) exp)
    {
    case 0:
        return VarDbl(1);
    case 1:
        return *this;
    default:
        break;
    }
    if ((exp > 0) && (std::ceil(exp) == std::floor(exp))) {
        std::vector<VarDbl> sCoeff(exp + 1);
        sCoeff[(size_t) exp] = 1;
        return polynominal(sCoeff);
    }
    std::vector<VarDbl> sTaylor;
    sTaylor.reserve(MAX_ORDER_FOR_TAYLOR);
    sTaylor.push_back(VarDbl(std::pow(value(), exp)));
    sTaylor.push_back(VarDbl(exp));
    for (size_t i = 2; i < MAX_ORDER_FOR_TAYLOR; ++i) {
        sTaylor.push_back(sTaylor[i - 1] * ((exp + 1)/i - 1));
        if (sTaylor.back().value() == 0)
            break;
    }
    std::ostringstream os;
    os << "pow(" << *this << ", " << exp << ")";
    return taylor1d(os.str(), sTaylor, true, true);
}

} // namespace var_dbl
#endif // __ValDbl_h__
