/*
Best-effort double-precision interval arithmetic with ±ULP outward padding
after each operation. Not IEEE-1788 rigorous (no directed-rounding mode),
but tight enough for diagnostic comparison against variance arithmetic.

Conversion from VarDbl: [value - max(unc, ulp(value)), value + max(unc, ulp(value))].
Useful for wrapping IndexSin Quart twiddles (which store value with unc=0) into
proper intervals bounding the underlying double-precision sin/cos error.
*/
#ifndef __Interval_h__
#define __Interval_h__

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>

#include "ulp.h"
#include "VarDbl.h"

namespace var_dbl
{

class Interval {
    double _lo;
    double _hi;

    static double _safeUlp(double x) {
        const double u = ulp(x);
        return (u > 0.0) ? u : std::numeric_limits<double>::min();
    }
    static Interval _padded(double lo, double hi) {
        Interval out;
        out._lo = lo - _safeUlp(lo);
        out._hi = hi + _safeUlp(hi);
        return out;
    }

public:
    Interval() : _lo(0.0), _hi(0.0) {}
    Interval(double point) : _lo(point), _hi(point) {}
    Interval(double lo, double hi)
        : _lo(std::min(lo, hi)), _hi(std::max(lo, hi)) {}
    // From a VarDbl: half-width = max(uncertainty, ulp(value)) so degenerate
    // (unc=0) sources like IndexSin Quart still get a proper interval.
    explicit Interval(const VarDbl& v) {
        const double w = std::max(v.uncertainty(), _safeUlp(v.value()));
        _lo = v.value() - w;
        _hi = v.value() + w;
    }
    static Interval centered(double value, double unc) {
        if (unc < 0)
            throw std::invalid_argument("Interval::centered: unc must be >= 0");
        const double w = std::max(unc, _safeUlp(value));
        return Interval(value - w, value + w);
    }

    double lo()  const { return _lo; }
    double hi()  const { return _hi; }
    double mid() const { return 0.5 * (_lo + _hi); }
    double rad() const { return 0.5 * (_hi - _lo); }
    bool contains(double x) const { return _lo <= x && x <= _hi; }
    bool encloses(const Interval& o) const { return _lo <= o._lo && o._hi <= _hi; }

    Interval operator-() const { return Interval(-_hi, -_lo); }

    Interval& operator+=(const Interval& r) {
        *this = _padded(_lo + r._lo, _hi + r._hi);
        return *this;
    }
    Interval& operator-=(const Interval& r) {
        *this = _padded(_lo - r._hi, _hi - r._lo);
        return *this;
    }
    Interval& operator*=(const Interval& r) {
        const double a = _lo * r._lo, b = _lo * r._hi;
        const double c = _hi * r._lo, d = _hi * r._hi;
        *this = _padded(std::min(std::min(a, b), std::min(c, d)),
                        std::max(std::max(a, b), std::max(c, d)));
        return *this;
    }
    Interval& operator*=(double s) {
        const double a = _lo * s, b = _hi * s;
        *this = _padded(std::min(a, b), std::max(a, b));
        return *this;
    }
};

inline Interval operator+(Interval a, const Interval& b) { a += b; return a; }
inline Interval operator-(Interval a, const Interval& b) { a -= b; return a; }
inline Interval operator*(Interval a, const Interval& b) { a *= b; return a; }
inline Interval operator*(Interval a, double s)          { a *= s; return a; }
inline Interval operator*(double s, Interval a)          { a *= s; return a; }

inline std::ostream& operator<<(std::ostream& os, const Interval& iv) {
    return os << "[" << iv.lo() << ", " << iv.hi() << "]";
}

} // namespace var_dbl
#endif
