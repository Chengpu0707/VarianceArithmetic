/*
ImPrecise<T> is a template class implementing variance arithmetic in C++.
T is constrained to a floating-point type and is used for both the value
and the standard uncertainty (sqrt of variance) of each instance. VarDbl is
provided as a typedef for ImPrecise<double> to preserve the existing API.

VarDbl.h is header-only and needs no build system. It relies only on ulp.h
and Test.h (for assertEqual()), both also in the var_dbl namespace.
Division (operator/=) is defined in Taylor.h, as the expansion of 1/x.
*/
#ifndef __VarDbl_h__
#define __VarDbl_h__

#include <algorithm>
#include <cmath>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#if __cplusplus >= 201103L
#include <type_traits>
#endif
#if __cplusplus >= 202002L
#include <concepts>
#endif

#include "Test.h"
#include "ulp.h"

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

    explicit InitException(double value, double uncertainty, const std::string& what) :
        runtime_error(what), value(value), uncertainty(uncertainty)
    {}
};


IMP_CLASS_TMPL
class ImPrecise {
#if __cplusplus >= 201103L
    static_assert(std::is_floating_point<T>::value,
                  "ImPrecise<T>: T must be a floating-point type");
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

    // Store (value, |uncertainty|), or throw InitException if either is not finite.
    void init(T value, T uncertainty, const char* what);
    // Throw InitException with the message "<prefix><value>[, <uncertainty>])".
    // Building the message only on failure keeps the constructors cheap.
    static void throwInit(const char* prefix, T value, T uncertainty, bool withUncertainty);

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
    ImPrecise();
    ImPrecise(T value, T uncertainty);

    // conversion constructors
#if __cplusplus >= 202002L
    template <typename U> requires std::floating_point<U> ImPrecise(U value);
    template <typename U> requires std::integral<U> ImPrecise(U value);
#elif __cplusplus >= 201103L
    template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type = 0>
    ImPrecise(U value);
    template <typename U, typename std::enable_if<std::is_integral<U>::value, int>::type = 0>
    ImPrecise(U value);
#else
    template <typename U>
    ImPrecise(U value, typename _enable_if_c03<_is_floating_point_c03<U>::value, _C03FloatTag>::type = _C03FloatTag());
    template <typename U>
    ImPrecise(U value, typename _enable_if_c03<_is_integral_c03<U>::value, _C03IntTag>::type = _C03IntTag());
#endif

    // i/o
    std::string to_string() const;

    // Inline-friend i/o (single-instantiation, avoids template-friend pitfalls).
    friend std::ostream& operator<<(std::ostream& out, const ImPrecise& v) {
        out << std::scientific << v._value << '~' << v._uncertainty;
        return out;
    }
    friend std::istream& operator>>(std::istream& in, ImPrecise& v) {
        T value = T(0), uncertainty = T(0);
        char sep = '\0';
        if (!(in >> value >> sep >> uncertainty))
            return in;
        if (sep != '~') {
            std::ostringstream ss;
            ss << "Invalid separator for " << value << sep << uncertainty;
            throw std::invalid_argument(ss.str());
        }
        v.init(value, uncertainty, "operator>>");
        return in;
    }

    // +
    ImPrecise operator+(const ImPrecise& other) const;
    ImPrecise& operator+=(const ImPrecise& other);
    IMP_SCALAR_TMPL ImPrecise operator+(const U& other) const;
    IMP_SCALAR_TMPL ImPrecise& operator+=(const U& other);

    // -
    void negate() { _value = -_value; }
    ImPrecise operator-() const;
    ImPrecise operator-(const ImPrecise& other) const;
    ImPrecise& operator-=(const ImPrecise& other);
    IMP_SCALAR_TMPL ImPrecise operator-(const U& other) const;
    IMP_SCALAR_TMPL ImPrecise& operator-=(const U& other);

    // *
    ImPrecise operator*(const ImPrecise& other) const;
    ImPrecise& operator*=(const ImPrecise& other);
    IMP_SCALAR_TMPL ImPrecise operator*(const U& other) const;
    IMP_SCALAR_TMPL ImPrecise& operator*=(const U& other);

    // /
    ImPrecise operator/(const ImPrecise& other) const;
    ImPrecise& operator/=(const ImPrecise& other);
    IMP_SCALAR_TMPL ImPrecise operator/(const U& other) const;
    IMP_SCALAR_TMPL ImPrecise& operator/=(const U& other);

    // Symmetric inline-friend scalar operators (scalar on the left).
    IMP_SCALAR_TMPL friend ImPrecise operator+(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) + second;
    }
    IMP_SCALAR_TMPL friend ImPrecise operator-(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) - second;
    }
    IMP_SCALAR_TMPL friend ImPrecise operator*(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) * second;
    }
    IMP_SCALAR_TMPL friend ImPrecise operator/(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) / second;
    }

    // compare
    bool operator==(const ImPrecise& other) const;
    bool operator!=(const ImPrecise& other) const;
    bool operator<(const ImPrecise& other) const;
    bool operator>(const ImPrecise& other) const;
    bool operator<=(const ImPrecise& other) const;
    bool operator>=(const ImPrecise& other) const;
    IMP_SCALAR_TMPL bool operator==(const U& other) const;
    IMP_SCALAR_TMPL bool operator!=(const U& other) const;
    IMP_SCALAR_TMPL bool operator<(const U& other) const;
    IMP_SCALAR_TMPL bool operator>(const U& other) const;
    IMP_SCALAR_TMPL bool operator<=(const U& other) const;
    IMP_SCALAR_TMPL bool operator>=(const U& other) const;

    IMP_SCALAR_TMPL friend bool operator==(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) == second;
    }
    IMP_SCALAR_TMPL friend bool operator!=(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) != second;
    }
    IMP_SCALAR_TMPL friend bool operator<(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) < second;
    }
    IMP_SCALAR_TMPL friend bool operator>(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) > second;
    }
    IMP_SCALAR_TMPL friend bool operator<=(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) <= second;
    }
    IMP_SCALAR_TMPL friend bool operator>=(const U& first, const ImPrecise& second) {
        return ImPrecise(static_cast<T>(first)) >= second;
    }

    static void assertEqual(const ImPrecise& var, T value, T uncertainty,
                            const std::string& msg = "",
                            T valueDelta = T(0), T uncertaintyDelta = T(0));

};


// ---------- Inline implementations ----------

IMP_CLASS_TMPL
inline void ImPrecise<T>::init(T value, T uncertainty, const char* what) {
    if (!std::isfinite(value) || !std::isfinite(uncertainty))
        throw InitException(static_cast<double>(value), static_cast<double>(uncertainty), what);
    _value = value;
    _uncertainty = std::abs(uncertainty);
}

IMP_CLASS_TMPL
inline void ImPrecise<T>::throwInit(const char* prefix, T value, T uncertainty, bool withUncertainty) {
    std::ostringstream ss;
    ss << prefix << value;
    if (withUncertainty)
        ss << ", " << uncertainty;
    ss << ")";
    throw InitException(static_cast<double>(value), static_cast<double>(uncertainty), ss.str());
}


IMP_CLASS_TMPL
inline ImPrecise<T>::ImPrecise() : _value(T(0)), _uncertainty(T(0)) {
}

IMP_CLASS_TMPL
inline ImPrecise<T>::ImPrecise(T value, T uncertainty)
{
    if (!std::isfinite(value) || !std::isfinite(uncertainty))
        throwInit("ImPrecise( ", value, uncertainty, true);
    _value = value;
    _uncertainty = std::abs(uncertainty);
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U>
inline T ImPrecise<T>::ulp(U value) {
#elif __cplusplus >= 201103L
template <typename T>
template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type>
inline T ImPrecise<T>::ulp(U value) {
#else
template <typename T>
template <typename U>
inline typename _enable_if_c03<_is_floating_point_c03<U>::value, T>::type
ImPrecise<T>::ulp(U value) {
#endif
    return static_cast<T>(var_dbl::ulp(value, PRECISE_SIGNIFICAND_TAIL_BITS_VALUE))
            * DEVIATION_OF_LSB();
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U>
inline ImPrecise<T>::ImPrecise(U value) {
#elif __cplusplus >= 201103L
template <typename T>
template <typename U, typename std::enable_if<std::is_floating_point<U>::value, int>::type>
inline ImPrecise<T>::ImPrecise(U value) {
#else
template <typename T>
template <typename U>
inline ImPrecise<T>::ImPrecise(U value, typename _enable_if_c03<_is_floating_point_c03<U>::value, _C03FloatTag>::type) {
#endif
    const T v = static_cast<T>(value);
    if (!std::isfinite(v))
        throwInit("ImPrecise(fp ", v, T(0), false);
    _value = v;
    _uncertainty = std::abs(ImPrecise<T>::ulp(v));
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::integral<U>
inline ImPrecise<T>::ImPrecise(U value) {
#elif __cplusplus >= 201103L
template <typename T>
template <typename U, typename std::enable_if<std::is_integral<U>::value, int>::type>
inline ImPrecise<T>::ImPrecise(U value) {
#else
template <typename T>
template <typename U>
inline ImPrecise<T>::ImPrecise(U value, typename _enable_if_c03<_is_integral_c03<U>::value, _C03IntTag>::type) {
#endif
    const T v = static_cast<T>(value);
    if (!std::isfinite(v))
        throwInit("ImPrecise(i ", v, T(0), false);
    _value = v;
    _uncertainty = static_cast<T>(std::abs(var_dbl::ulp(value)));
}


IMP_CLASS_TMPL
inline std::string ImPrecise<T>::to_string() const {
    std::ostringstream os;
    os << std::scientific << value() << '~' << uncertainty();
    return os.str();
}


IMP_CLASS_TMPL
inline ImPrecise<T> ImPrecise<T>::operator-() const {
    ImPrecise ret(*this);
    ret.negate();
    return ret;
}

IMP_CLASS_TMPL
inline ImPrecise<T>& ImPrecise<T>::operator+=(const ImPrecise& other) {
    const T u = (this->_uncertainty == T(0)) ? other._uncertainty :
                (other._uncertainty == T(0)) ? this->_uncertainty :
                std::sqrt(variance() + other.variance());
    init(_value + other._value, u, "+=");
    return *this;
}

IMP_CLASS_TMPL
inline ImPrecise<T>& ImPrecise<T>::operator-=(const ImPrecise& other) {
    const T u = (this->_uncertainty == T(0)) ? other._uncertainty :
                (other._uncertainty == T(0)) ? this->_uncertainty :
                std::sqrt(variance() + other.variance());
    init(_value - other._value, u, "-=");
    return *this;
}

IMP_CLASS_TMPL
inline ImPrecise<T>& ImPrecise<T>::operator*=(const ImPrecise& other) {
    const T u = std::sqrt(
            this->variance() * other.value() * other.value() +
            other.variance() * value() * value() +
            this->variance() * other.variance());
    const T product = value() * other.value();
    if (u <= T(0)) {
        // Both operands are precise. If they are integers whose exact product
        // exceeds the significand of T, the rounding of the product is its only
        // uncertainty, so multiply exactly as integers. The checks keep every
        // cast to long long exact and the integer product from overflowing.
        const T maxSig = static_cast<T>(ImPrecise<T>::MAX_SIGNIFICAND());
        const T a = value(), b = other.value();
        if (std::abs(a) <= maxSig && std::floor(a) == a &&
            std::abs(b) <= maxSig && std::floor(b) == b &&
            maxSig < std::abs(product) &&
            std::abs(product) < static_cast<T>(std::numeric_limits<long long>::max())) {
            *this = ImPrecise<T>(static_cast<long long>(a) * static_cast<long long>(b));
            return *this;
        }
    }
    init(product, u, "*=");
    return *this;
}

// ImPrecise<T>& ImPrecise<T>::operator/=(const ImPrecise& other) in Taylor.h


IMP_CLASS_TMPL
inline ImPrecise<T> ImPrecise<T>::operator+(const ImPrecise& other) const {
    ImPrecise res(*this);
    res += other;
    return res;
}

IMP_CLASS_TMPL
inline ImPrecise<T> ImPrecise<T>::operator-(const ImPrecise& other) const {
    ImPrecise res(*this);
    res -= other;
    return res;
}

IMP_CLASS_TMPL
inline ImPrecise<T> ImPrecise<T>::operator*(const ImPrecise& other) const {
    ImPrecise res(*this);
    res *= other;
    return res;
}

IMP_CLASS_TMPL
inline ImPrecise<T> ImPrecise<T>::operator/(const ImPrecise& other) const {
    ImPrecise res(*this);
    res /= other;
    return res;
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T> ImPrecise<T>::operator+(const U& other) const {
#else
template <typename T>
template <typename U>
inline ImPrecise<T> ImPrecise<T>::operator+(const U& other) const {
#endif
    return *this + ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T> ImPrecise<T>::operator-(const U& other) const {
#else
template <typename T>
template <typename U>
inline ImPrecise<T> ImPrecise<T>::operator-(const U& other) const {
#endif
    return *this - ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T> ImPrecise<T>::operator*(const U& other) const {
#else
template <typename T>
template <typename U>
inline ImPrecise<T> ImPrecise<T>::operator*(const U& other) const {
#endif
    return *this * ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T> ImPrecise<T>::operator/(const U& other) const {
#else
template <typename T>
template <typename U>
inline ImPrecise<T> ImPrecise<T>::operator/(const U& other) const {
#endif
    return *this / ImPrecise<T>(other);
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T>& ImPrecise<T>::operator+=(const U& other) {
#else
template <typename T>
template <typename U>
inline ImPrecise<T>& ImPrecise<T>::operator+=(const U& other) {
#endif
    return *this += ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T>& ImPrecise<T>::operator-=(const U& other) {
#else
template <typename T>
template <typename U>
inline ImPrecise<T>& ImPrecise<T>::operator-=(const U& other) {
#endif
    return *this -= ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T>& ImPrecise<T>::operator*=(const U& other) {
#else
template <typename T>
template <typename U>
inline ImPrecise<T>& ImPrecise<T>::operator*=(const U& other) {
#endif
    return *this *= ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline ImPrecise<T>& ImPrecise<T>::operator/=(const U& other) {
#else
template <typename T>
template <typename U>
inline ImPrecise<T>& ImPrecise<T>::operator/=(const U& other) {
#endif
    return *this /= ImPrecise<T>(other);
}


IMP_CLASS_TMPL
inline bool ImPrecise<T>::operator==(const ImPrecise& other) const {
    ImPrecise res(*this);
    res -= other;
    return std::abs(res.value()) <= (ImPrecise<T>::BINDING_FOR_EQUAL() * res.uncertainty());
}

IMP_CLASS_TMPL
inline bool ImPrecise<T>::operator!=(const ImPrecise& other) const {
    return !(*this == other);
}

IMP_CLASS_TMPL
inline bool ImPrecise<T>::operator<(const ImPrecise& other) const {
    if (*this == other)
        return false;
    return this->value() < other.value();
}

IMP_CLASS_TMPL
inline bool ImPrecise<T>::operator>(const ImPrecise& other) const {
    if (*this == other)
        return false;
    return this->value() > other.value();
}

IMP_CLASS_TMPL
inline bool ImPrecise<T>::operator<=(const ImPrecise& other) const {
    if (*this == other)
        return true;
    return this->value() < other.value();
}

IMP_CLASS_TMPL
inline bool ImPrecise<T>::operator>=(const ImPrecise& other) const {
    if (*this == other)
        return true;
    return this->value() > other.value();
}


#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool ImPrecise<T>::operator==(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool ImPrecise<T>::operator==(const U& other) const {
#endif
    return *this == ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool ImPrecise<T>::operator!=(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool ImPrecise<T>::operator!=(const U& other) const {
#endif
    return *this != ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool ImPrecise<T>::operator<(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool ImPrecise<T>::operator<(const U& other) const {
#endif
    return *this < ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool ImPrecise<T>::operator>(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool ImPrecise<T>::operator>(const U& other) const {
#endif
    return *this > ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool ImPrecise<T>::operator<=(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool ImPrecise<T>::operator<=(const U& other) const {
#endif
    return *this <= ImPrecise<T>(other);
}

#if __cplusplus >= 202002L
template <typename T> requires std::floating_point<T>
template <typename U> requires std::floating_point<U> || std::integral<U>
inline bool ImPrecise<T>::operator>=(const U& other) const {
#else
template <typename T>
template <typename U>
inline bool ImPrecise<T>::operator>=(const U& other) const {
#endif
    return *this >= ImPrecise<T>(other);
}


IMP_CLASS_TMPL
inline void ImPrecise<T>::assertEqual(const ImPrecise& var, T value, T uncertainty,
            const std::string& msg,
            T deltaValue, T deltaUncertainty) {
    std::ostringstream oss;
    if (!msg.empty())
        oss << msg << ": ";
    oss << var.value() << "~" << var.uncertainty() << " vs " << value << "~" << uncertainty;
    if (deltaValue == T(0))
        deltaValue = std::max(ImPrecise<T>::ulp(var.value()), ImPrecise<T>::ulp(value));
    test::assertAlmostEqual(static_cast<double>(var.value()),
                            static_cast<double>(value),
                            static_cast<double>(deltaValue), oss.str());
    if (deltaUncertainty == T(0))
        deltaUncertainty = std::max(ImPrecise<T>::ulp(var.uncertainty()), ImPrecise<T>::ulp(uncertainty));
    test::assertAlmostEqual(static_cast<double>(var.uncertainty()),
                            static_cast<double>(uncertainty),
                            static_cast<double>(deltaUncertainty), oss.str());
}


// Preserve the existing API: VarDbl is ImPrecise<double>.
typedef ImPrecise<double> VarDbl;


} // namespace var_dbl
#endif // __ValDbl_h__
