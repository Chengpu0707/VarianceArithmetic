#include <cmath>
#include <exception>
#include <format>   // not supported in gcc 13.2.0
#include <sstream>
#include <string>

#include "Test.h"

#ifndef __ValDbl_h__
#define __ValDbl_h__

class ValueError : public std::exception {
public:
    const double value;
    const std::string what;

    explicit ValueError(double value, const std::string what) : 
        value(value), what(what)
    {
    }
};


class UncertaintyError : public std::exception {
public:
    const double value;
    const double uncertainty;
    const std::string what;

    explicit UncertaintyError(double value, double uncertainty, const std::string what) : 
        value(value), uncertainty(uncertainty), what(what)
    {
    }
};

class VarDbl {
private:
    double _value = 0;
    double _variance = 0;

    // assume uniform distribution within ulp()
    constexpr static const double DEVIATION_OF_LSB = 1.0 / sqrt(3);

public:
    double value() const { return _value; }
    double uncertainty() const { return sqrt(_variance); }

    // constructors
    explicit VarDbl(double value, double uncertainty);
        // uncertainty is limited between sqrt(std::numeric_limits<double>::mim()) and sqrt(std::numeric_limits<double>::max())
    explicit VarDbl(double value);
        // assume ulp as uncertainty
    explicit VarDbl(long long value) noexcept;
        // may loss resolution
    explicit VarDbl(long value) : VarDbl((long long) value) {}
    explicit VarDbl(int value) : VarDbl((long long) value) {}

private:
    void init(double value, double uncertainty, const std::string what) 
    {
        if (!std::isfinite(value))
            throw ValueError(value, what);
        const double variance = uncertainty * uncertainty;
        if (!std::isfinite(variance))
            throw UncertaintyError(value, uncertainty, what);
        _value = value;
        _variance = variance;
    }
};


inline VarDbl::VarDbl(double value, double uncertainty) 
{
    std::ostringstream ss;
    ss << "VarDbl(double " << value << ", double " << uncertainty << ")";
    init(value, uncertainty, ss.str());
}

inline VarDbl::VarDbl(double value) 
{
    std::ostringstream ss;
    ss << "VarDbl(double " << value << ")";
    init(value, Test::ulp(value)*DEVIATION_OF_LSB, ss.str());
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
    init(value, Test::ulp(value)*DEVIATION_OF_LSB, ss.str());
}


#endif // __ValDbl_h__