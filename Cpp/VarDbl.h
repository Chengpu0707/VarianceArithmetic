/*
Imprecise<T> is a template class implementing variance arithmetic in C++.
T is constrained to a floating-point type and is used for both the value
and the standard uncertainty (sqrt of variance) of each instance. VarDbl is
provided as a typedef for Imprecise<double> to preserve the existing API.

VarDbl.h is intended to be fully inline and deplored along without any
making system. It relies only on Moment.h, which is also in the var_dbl
namespace.
*/
#include <cmath>
#include <limits>

#include "Test.h"
#include "ulp.h"

#ifndef __ValDbl_h__
#define __ValDbl_h__

// Class template constraint: T must be a floating-point type.
#if __cplusplus >= 202002L
#define IMP_CLASS_TMPL template <typename T> requires std::floating_point<T>
#else
#define IMP_CLASS_TMPL template <typename T>
#endif

// Scalar operand constraint inside member templates: floating-point OR integral.
// Uses parameter name U to avoid shadowing the class template's T.
#if __cplusplus >= 202002L
#define IMP_SCALAR_TMPL template <typename U> requires std::floating_point<U> || std::integral<U>
#else
#define IMP_SCALAR_TMPL template <typename U>
#endif

namespace var_dbl
{


struct InitException : public std::runtime_error
{
    const double value;
    const double uncertainty;

    explicit InitException(double value, double uncertainty, const std::string what) :
        runtime_error(what), value(value), uncertainty(uncertainty)
    {}
};


IMP_CLASS_TMPL
class Imprecise {
#if __cplusplus >= 201103L
    static_assert(std::is_floating_point<T>::value,
                  "Imprecise<T>: T must be a floating-point type");
#endif

    // If the last 23 bits of significand are zero, the value is 2's fractional.
    static const long PRECISE_SIGNIFICAND_TAIL_BITS_VALUE = 23;

    // Max integer significand representable in T (= 2^digits - 1).
    static long long MAX_SIGNIFICAND() {
        return (1LL << std::numeric_limits<T>::digits) - 1LL;
    }

    // Deviation factor for a uniform distribution within ULP: 1/sqrt(3).
    static T DEVIATION_OF_LSB() {
        return T(1) / std::sqrt(T(3));
    }

    // z-score for 50% probability of equal.
    static T BINDING_FOR_EQUAL() {
        return T(0.67448975);
    }

#if __cplusplus >= 201103L
    T _value = T(0);
    T _uncertainty = T(0);
#else
    T _value;
    T _uncertainty;
#endif

    void init(T value, T uncertainty, const std::string what);

public:
    // ulp helper — returns T-typed standard uncertainty for value.
#if __cplusplus >= 202002L
    template <typename U> requires std::floating_point<U> static T ulp(U value);
#elif __cplusplus >= 201103L
    template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type = 0>
    static T ulp(U value);
#else
    template <typename U>
    static typename _enable_if_c03<_is_floating_point_c03<U>::value, T>::type ulp(U value);
#endif

    T value() const { return _value; }
    T uncertainty() const { return _uncertainty; }
    T variance() const { return _uncertainty * _uncertainty; }

    // constructors
    Imprecise();
    Imprecise(const Imprecise& other);
    Imprecise(T value, T uncertainty);

    // conversion constructors
#if __cplusplus >= 202002L
    template <typename U> requires std::floating_point<U> Imprecise(U value);
    template <typename U> requires std::integral<U> Imprecise(U value);
#elif __cplusplus >= 201103L
    template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type = 0>
    Imprecise(U value);
    template <typename U, typename std::enable_if<std::is_integral<U>::value, int>::type = 0>
    Imprecise(U value);
#else
    template <typename U>
    Imprecise(U value, typename _enable_if_c03<_is_floating_point_c03<U>::value, _C03FloatTag>::type = _C03FloatTag());
    template <typename U>
    Imprecise(U value, typename _enable_if_c03<_is_integral_c03<U>::value, _C03IntTag>::type = _C03IntTag());
#endif

    // i/o
    std::string to_string() const;

    // Inline-friend i/o (single-instantiation, avoids template-friend pitfalls).
    friend std::ostream& operator<<(std::ostream& out, const Imprecise& v) {
        out << std::scientific << v._value << '~' << v._uncertainty;
        return out;
    }
    friend std::istream& operator>>(std::istream& in, Imprecise& v) {
        T value, uncertainty;
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

    // +
    Imprecise operator+(const Imprecise& other) const;
    Imprecise operator+=(const Imprecise& other);
    IMP_SCALAR_TMPL Imprecise operator+(const U& other) const;
    IMP_SCALAR_TMPL Imprecise operator+=(const U& other);

    // -
    void negate() { _value = -_value; }
    Imprecise operator-() const;
    Imprecise operator-(const Imprecise& other) const;
    Imprecise operator-=(const Imprecise& other);
    IMP_SCALAR_TMPL Imprecise operator-(const U& other) const;
    IMP_SCALAR_TMPL Imprecise operator-=(const U& other);

    // *
    Imprecise operator*(const Imprecise& other) const;
    Imprecise operator*=(const Imprecise& other);
    IMP_SCALAR_TMPL Imprecise operator*(const U& other) const;
    IMP_SCALAR_TMPL Imprecise operator*=(const U& other);

    // /
    Imprecise operator/(const Imprecise& other) const;
    Imprecise operator/=(const Imprecise& other);
    IMP_SCALAR_TMPL Imprecise operator/(const U& other) const;
    IMP_SCALAR_TMPL Imprecise operator/=(const U& other);

    // Symmetric inline-friend scalar operators (scalar on the left).
    IMP_SCALAR_TMPL friend Imprecise operator+(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) + second;
    }
    IMP_SCALAR_TMPL friend Imprecise operator-(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) - second;
    }
    IMP_SCALAR_TMPL friend Imprecise operator*(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) * second;
    }
    IMP_SCALAR_TMPL friend Imprecise operator/(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) / second;
    }

    // compare
    bool operator==(const Imprecise& other) const;
    bool operator!=(const Imprecise& other) const;
    bool operator<(const Imprecise& other) const;
    bool operator>(const Imprecise& other) const;
    bool operator<=(const Imprecise& other) const;
    bool operator>=(const Imprecise& other) const;
    IMP_SCALAR_TMPL bool operator==(const U& other) const;
    IMP_SCALAR_TMPL bool operator!=(const U& other) const;
    IMP_SCALAR_TMPL bool operator<(const U& other) const;
    IMP_SCALAR_TMPL bool operator>(const U& other) const;
    IMP_SCALAR_TMPL bool operator<=(const U& other) const;
    IMP_SCALAR_TMPL bool operator>=(const U& other) const;

    IMP_SCALAR_TMPL friend bool operator==(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) == second;
    }
    IMP_SCALAR_TMPL friend bool operator!=(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) != second;
    }
    IMP_SCALAR_TMPL friend bool operator<(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) < second;
    }
    IMP_SCALAR_TMPL friend bool operator>(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) > second;
    }
    IMP_SCALAR_TMPL friend bool operator<=(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) <= second;
    }
    IMP_SCALAR_TMPL friend bool operator>=(const U& first, const Imprecise& second) {
        return Imprecise(static_cast<T>(first)) >= second;
    }

    static void assertEqual(const Imprecise& var, T value, T uncertainty,
                            const std::string& msg = "",
                            T valueDelta = T(0), T uncertaintyDelta = T(0));

};


// ---------- Inline implementations ----------

IMP_CLASS_TMPL
inline void Imprecise<T>::init(T value, T uncertainty, const std::string what) {
    if (!std::isfinite(value) || !std::isfinite(uncertainty))
        throw InitException(static_cast<double>(value), static_cast<double>(uncertainty), what);
    _value = value;
    _uncertainty = std::abs(uncertainty);
}


IMP_CLASS_TMPL
inline Imprecise<T>::Imprecise() {
    _value = T(0);
    _uncertainty = T(0);
}

IMP_CLASS_TMPL
inline Imprecise<T>::Imprecise(const Imprecise& other) {
    _value = other._value;
    _uncertainty = other._uncertainty;
}

IMP_CLASS_TMPL
inline Imprecise<T>::Imprecise(T value, T uncertainty)
{
    std::ostringstream ss;
    ss << "Imprecise( " << value << ", " << uncertainty << ")";
    init(value, uncertainty, ss.str());
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U>
inline T Imprecise<T>::ulp(U value) {
#elif __cplusplus >= 201103L
template <typename T>
template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type>
inline T Imprecise<T>::ulp(U value) {
#else
template <typename T>
template <typename U>
inline typename _enable_if_c03<_is_floating_point_c03<U>::value, T>::type
Imprecise<T>::ulp(U value) {
#endif
    return static_cast<T>(var_dbl::ulp(value, PRECISE_SIGNIFICAND_TAIL_BITS_VALUE))
            * DEVIATION_OF_LSB();
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U>
inline Imprecise<T>::Imprecise(U value) {
#elif __cplusplus >= 201103L
template <typename T>
template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type>
inline Imprecise<T>::Imprecise(U value) {
#else
template <typename T>
template <typename U>
inline Imprecise<T>::Imprecise(U value, typename _enable_if_c03<_is_floating_point_c03<U>::value, _C03FloatTag>::type) {
#endif
    std::ostringstream ss;
    ss << "Imprecise(fp " << value << ")";
    const T v = static_cast<T>(value);
    const T u = Imprecise<T>::ulp(v);
    init(v, u, ss.str());
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::integral<U>
inline Imprecise<T>::Imprecise(U value) {
#elif __cplusplus >= 201103L
template <typename T>
template <typename U, typename std::enable_if<std::is_integral<U>::value, int>::type>
inline Imprecise<T>::Imprecise(U value) {
#else
template <typename T>
template <typename U>
inline Imprecise<T>::Imprecise(U value, typename _enable_if_c03<_is_integral_c03<U>::value, _C03IntTag>::type) {
#endif
    std::ostringstream ss;
    ss << "Imprecise(i " << value << ")";
    const T v = static_cast<T>(value);
    const T u = static_cast<T>(std::abs(var_dbl::ulp(value)));
    init(v, u, ss.str());
}


IMP_CLASS_TMPL
inline std::string Imprecise<T>::to_string() const {
    std::ostringstream os;
    os << std::scientific << value() << '~' << uncertainty();
    return os.str();
}


IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator-() const {
    Imprecise ret(*this);
    ret.negate();
    return ret;
}

IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator+=(const Imprecise& other) {
    const T u = (this->_uncertainty == T(0)) ? other._uncertainty :
                (other._uncertainty == T(0)) ? this->_uncertainty :
                std::sqrt(variance() + other.variance());
    init(_value + other._value, u, "+=");
    return *this;
}

IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator-=(const Imprecise& other) {
    const T u = (this->_uncertainty == T(0)) ? other._uncertainty :
                (other._uncertainty == T(0)) ? this->_uncertainty :
                std::sqrt(variance() + other.variance());
    init(_value - other._value, u, "-=");
    return *this;
}

IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator*=(const Imprecise& other) {
    const T u = std::sqrt(
            this->variance() * other.value() * other.value() +
            other.variance() * value() * value() +
            this->variance() * other.variance());
    if (u <= T(0)) {
        const long long val = static_cast<long long>(value()) * static_cast<long long>(other.value());
        const long long maxSig = Imprecise<T>::MAX_SIGNIFICAND();
        if (maxSig < std::abs(val) &&
            value() < static_cast<T>(maxSig) && other.value() < static_cast<T>(maxSig)) {
            const Imprecise<T> v(val);
            this->_value = v._value;
            this->_uncertainty = v._uncertainty;
        } else
            init(value() * other.value(), u, "*=");
    } else
        init(value() * other.value(), u, "*=");
    return *this;
}

// Imprecise<T> Imprecise<T>::operator/=(const Imprecise& other) in Taylor.h


IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator+(const Imprecise& other) const {
    Imprecise res(*this);
    res += other;
    return res;
}

IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator-(const Imprecise& other) const {
    Imprecise res(*this);
    res -= other;
    return res;
}

IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator*(const Imprecise& other) const {
    Imprecise res(*this);
    res *= other;
    return res;
}

IMP_CLASS_TMPL
inline Imprecise<T> Imprecise<T>::operator/(const Imprecise& other) const {
    Imprecise res(*this);
    res /= other;
    return res;
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator+(const U& other) const {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator+(const U& other) const {
#endif
    return *this + Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator-(const U& other) const {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator-(const U& other) const {
#endif
    return *this - Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator*(const U& other) const {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator*(const U& other) const {
#endif
    return *this * Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator/(const U& other) const {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator/(const U& other) const {
#endif
    return *this / Imprecise<T>(other);
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator+=(const U& other) {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator+=(const U& other) {
#endif
    return *this += Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator-=(const U& other) {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator-=(const U& other) {
#endif
    return *this -= Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator*=(const U& other) {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator*=(const U& other) {
#endif
    return *this *= Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline Imprecise<T> Imprecise<T>::operator/=(const U& other) {
#else
template <typename T>
template <typename U>
inline Imprecise<T> Imprecise<T>::operator/=(const U& other) {
#endif
    return *this /= Imprecise<T>(other);
}


IMP_CLASS_TMPL
inline bool Imprecise<T>::operator==(const Imprecise& other) const {
    Imprecise res(*this);
    res -= other;
    return std::abs(res.value()) <= (Imprecise<T>::BINDING_FOR_EQUAL() * res.uncertainty());
}

IMP_CLASS_TMPL
inline bool Imprecise<T>::operator!=(const Imprecise& other) const {
    return !(*this == other);
}

IMP_CLASS_TMPL
inline bool Imprecise<T>::operator<(const Imprecise& other) const {
    if (*this == other)
        return false;
    return this->value() < other.value();
}

IMP_CLASS_TMPL
inline bool Imprecise<T>::operator>(const Imprecise& other) const {
    if (*this == other)
        return false;
    return this->value() > other.value();
}

IMP_CLASS_TMPL
inline bool Imprecise<T>::operator<=(const Imprecise& other) const {
    if (*this == other)
        return true;
    return this->value() < other.value();
}

IMP_CLASS_TMPL
inline bool Imprecise<T>::operator>=(const Imprecise& other) const {
    if (*this == other)
        return true;
    return this->value() > other.value();
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool Imprecise<T>::operator==(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool Imprecise<T>::operator==(const U& other) const {
#endif
    return *this == Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool Imprecise<T>::operator!=(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool Imprecise<T>::operator!=(const U& other) const {
#endif
    return *this != Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool Imprecise<T>::operator<(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool Imprecise<T>::operator<(const U& other) const {
#endif
    return *this < Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool Imprecise<T>::operator>(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool Imprecise<T>::operator>(const U& other) const {
#endif
    return *this > Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool Imprecise<T>::operator<=(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool Imprecise<T>::operator<=(const U& other) const {
#endif
    return *this <= Imprecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool Imprecise<T>::operator>=(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool Imprecise<T>::operator>=(const U& other) const {
#endif
    return *this >= Imprecise<T>(other);
}


IMP_CLASS_TMPL
inline void Imprecise<T>::assertEqual(const Imprecise& var, T value, T uncertainty,
            const std::string& msg,
            T deltaValue, T deltaUncertainty) {
    std::ostringstream oss;
    if (!msg.empty())
        oss << msg << ": ";
    oss << var.value() << "~" << var.uncertainty() << " vs " << value << "~" << uncertainty;
    if (deltaValue == T(0))
        deltaValue = std::max(Imprecise<T>::ulp(var.value()), Imprecise<T>::ulp(value));
    test::assertAlmostEqual(static_cast<double>(var.value()),
                            static_cast<double>(value),
                            static_cast<double>(deltaValue), oss.str());
    if (deltaUncertainty == T(0))
        deltaUncertainty = std::max(Imprecise<T>::ulp(var.uncertainty()), Imprecise<T>::ulp(uncertainty));
    test::assertAlmostEqual(static_cast<double>(var.uncertainty()),
                            static_cast<double>(uncertainty),
                            static_cast<double>(deltaUncertainty), oss.str());
}


// Preserve the existing API: VarDbl is Imprecise<double>.
typedef Imprecise<double> VarDbl;


} // namespace var_dbl
#endif // __ValDbl_h__
