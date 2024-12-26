/*
VarDbl is a class to implement variance arithmetic in C++

VarDbl.h is intended to be fully inline and deplored along without any making system.
It relies only on Momentum.h, which is also in the var_dbl namespace.
*/
#include <cmath>
#include <exception>
#include <functional>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "Momentum.h"
#include "test.h"
#include "ulp.h"

#ifndef __ValDbl_h__
#define __ValDbl_h__
namespace var_dbl 
{

struct VarDbl;

using UnionVector = std::variant<std::vector<VarDbl>, std::vector<double>>;


struct VarDbl { 
    constexpr static const double DEVIATION_OF_LSB = 1.0 / sqrt(3);
        // assume uniform distribution within ulp()
    constexpr static const double BINDING_FOR_EQUAL = 0.67448975;
        // z for 50% probability of equal

    // Taylor expansion
    constexpr static const int BINDING_FOR_TAYLOR = 5;
    constexpr static const double VARIANCE_THRESHOLD = 1.0 /BINDING_FOR_TAYLOR /BINDING_FOR_TAYLOR;
    constexpr static const int MIN_MONOTONIC_COUNT = 20;
    constexpr static const double TAU = 7.18e-7;  
    constexpr static const long PRECISE_SIGNIFICAND_TAIL_BITS = 23;

    constexpr static const char* INPUT_HEADER = 
            "name\tvalue\tuncertainty\tvariance\tinPrec\toutPrec\tBinding\tMaxOrder\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive";
    constexpr static const char* EXPANSION_HEADER = 
            "2n\tmonotonics\tExponent\tMomentum\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty";
    constexpr static const char* OUTPUT_HEADER =
            "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty";            

    template <typename T> requires std::floating_point<T> static double ulp(T value);

private:    
    static const NormalMomentum _momentum;
    double _value = 0;
    double _variance = 0;

    void init(double value, double variance, const std::string what);

    VarDbl taylor1d(std::string name, const UnionVector s1dTaylor, bool inPrec, bool outPrec, 
                    const char* const dumpPath=nullptr,
                    bool checkMonotonic=true, bool checkStability=true, bool checkReliablity=true, bool checkPositive=true
                   ) const;

    VarDbl polynominal(UnionVector sCoeff, const char* const dumpPath=nullptr) const;

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
    VarDbl exp(const char* const dumpPath=nullptr) const;
    VarDbl log(const char* const dumpPath=nullptr) const;
    VarDbl sin(const char* const dumpPath=nullptr) const;
    VarDbl pow(double exp, const char* const dumpPath=nullptr) const;


    /*
        1d Taylor expansion for differential series {s1dTaylor} at {input}, with {name} for logging.
        When {inPrec} is true, calculate Taylor expnasion against the precision of {input}.
        When {outPrec} is true, the result of the Taylor expnasion is the precision.
        s1dTaylor[n] should already normalized by /n!. 
        The max order of expansion is {maxOrder}, which should not exceed momentum.Normal.MAX_ORDER

        When {checkMonotonic} is true, raise {NotMonotonicException} if 
            after full expansion, the monotonic count is still less than {MIN_MONOTONIC_COUNT}.
        It should always be True.

        When {checkStability} is true, raise {NotStableException} if 
            after full expansion, the value for the last expansion term is more than 
            TAU-fold of the expansion uncertainty.
        It should always be True.

        When {checkReliablity} is true, raise {NotReliableException} if
            the precision of the result variance is more than 1/5, in which 5 is the binding factor.
        It should always be True.

        When {checkPositive} is true, raise {NotPosive} if
            the expansion variance at any order becomes negative
        It should always be True.

        Both the result value and variance are guaranteed to be finite, otherwise
            raise {NotFiniteException} 

        Dump the expansion to {dumpPath} when it is provided.
        {dumpPath} can be read back and tested using verifyDumpFile()
     */
    VarDbl taylor1d(std::string name, const std::vector<double>& s1dTaylor, bool inPrec, bool outPrec,
                    const char* const dumpPath=nullptr,
                    bool checkMonotonic=true, bool checkStability=true, bool checkReliablity=true, bool checkPositive=true
                   ) const;
    VarDbl taylor1d(std::string name, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
                    const char* const dumpPath=nullptr,
                    bool checkMonotonic=true, bool checkStability=true, bool checkReliablity=true, bool checkPositive=true
                   ) const;

    /*
     1d polynominal with coeeficient {sCoeff}
    */
    VarDbl polynominal(const std::vector<double>& sCoeff,
                       const char* const dumpPath=nullptr) const;
    VarDbl polynominal(const std::vector<VarDbl>& sCoeff,
                       const char* const dumpPath=nullptr) const;

};

const NormalMomentum VarDbl::_momentum(VarDbl::BINDING_FOR_TAYLOR);



struct InitException : public std::runtime_error
{
    const double value;
    const double variance;

    explicit InitException(double value, double variance, const std::string what) : 
        runtime_error(what), value(value), variance(variance)
    {}
};


struct TaylorIdException : public std::runtime_error 
{
    const VarDbl input;
    const bool inPrec, outPrec;
    const UnionVector s1dTaylor;
    const VarDbl value;
    const VarDbl variance;
    const size_t n;
    const size_t monotonics;
    const VarDbl newValue;
    const VarDbl newVariance;

    explicit TaylorIdException(const std::string what, 
                const VarDbl& input, const UnionVector& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance, 
                size_t n, size_t monotonics, const VarDbl& newValue, const VarDbl& newVariance) :
            runtime_error(what), input(input), 
            s1dTaylor(s1dTaylor), inPrec(inPrec), outPrec(outPrec),
            value(value), variance(variance), n(n), monotonics(monotonics),
            newValue(newValue), newVariance(newVariance) {
    }

    virtual std::string name() const { return "TaylorIdException"; }
};

struct NotFiniteException : public TaylorIdException
{
    using TaylorIdException::TaylorIdException;
    std::string name() const override { return "NotFiniteException"; }
};

struct NotReliableException : public TaylorIdException
{
    using TaylorIdException::TaylorIdException;
    std::string name() const override { return "NotReliableException"; }
};

struct NotMonotonicException : public TaylorIdException 
{
    using TaylorIdException::TaylorIdException;
    std::string name() const override { return "NotMonotonicException"; }
};

struct NotStableException : public TaylorIdException 
{
    using TaylorIdException::TaylorIdException;
    std::string name() const override { return "NotStableException"; }
};

struct NotPositiveException : public TaylorIdException 
{
    using TaylorIdException::TaylorIdException;
    std::string name() const override { return "NotPositiveException"; }
};






inline void VarDbl::init(double value, double variance, const std::string what) {
    if (!std::isfinite(value) || !std::isfinite(variance))
        throw InitException(value, variance, what);
    _value = value;
    _variance = variance;
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

inline VarDbl VarDbl::taylor1d(
    std::string name, const std::vector<double>& s1dTaylor, bool inPrec, bool outPrec,
    const char* const dumpPath,
    bool checkMonotonic, bool checkStability, bool checkReliablity, bool checkPositive
    ) const
{
    return taylor1d(name, UnionVector(std::move(s1dTaylor)), inPrec, outPrec,
                    dumpPath,
                    checkMonotonic, checkStability, checkReliablity, checkPositive);
}


inline VarDbl VarDbl::taylor1d(
    std::string name, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
    const char* const dumpPath,
    bool checkMonotonic, bool checkStability, bool checkReliablity, bool checkPositive
    ) const
{
    return taylor1d(name, UnionVector(std::move(s1dTaylor)), inPrec, outPrec,
                    dumpPath,
                    checkMonotonic, checkStability, checkReliablity, checkPositive);
}


inline VarDbl VarDbl::taylor1d(
    std::string name, const UnionVector s1dTaylor, bool inPrec, bool outPrec,
    const char* const dumpPath,
    bool checkMonotonic, bool checkStability, bool checkReliablity, bool checkPositive
    ) const
{
    const bool isDbl = std::holds_alternative<std::vector<double>>(s1dTaylor);
    const size_t size = isDbl? std::get<std::vector<double>>(s1dTaylor).size() : std::get<std::vector<VarDbl>>(s1dTaylor).size();
    const double* pDbl = isDbl? std::get<std::vector<double>>(s1dTaylor).data() : nullptr;
    const VarDbl* pVar = isDbl? nullptr : std::get<std::vector<VarDbl>>(s1dTaylor).data();

    if (isDbl) {
        for (int j = 0; j < size; ++j) {
            if (!std::isfinite(pDbl[j])) {
                std::ostringstream oss;
                oss << "Taylor " << pDbl[j] << " at " << j << "/" << size;
                throw std::invalid_argument(oss.str());
            }
        }
    } else {
        for (int j = 0; j < size; ++j) {
            if (!std::isfinite(pVar[j].value()) || !std::isfinite(pVar[j].variance())) {
                std::ostringstream oss;
                oss << "Taylor " << pVar[j] << " at " << j << "/" << size;
                throw std::invalid_argument(oss.str());
            }
        }
    }
     if (variance() == 0) {
        if (isDbl)
            return VarDbl(pDbl[0]);
        else
            return pVar[0];
    }
    std::ofstream ofs(dumpPath);
    if (ofs) {
        ofs << INPUT_HEADER <<"\n";
        ofs << name << '\t' << value() << '\t' << uncertainty() << '\t' << variance() << '\t' 
            << inPrec << '\t' << outPrec << '\t' 
            << BINDING_FOR_TAYLOR << '\t' << _momentum.size() << '\t' 
            << checkMonotonic << '\t' << checkStability << '\t' << checkReliablity << '\t' << checkPositive<< '\n';
        ofs << "Index";
        for (size_t i = 0; i < size; ++i)
            ofs << "\t" << i;
        ofs << "\nTaylor";
        if (isDbl) {
            for (size_t i = 0; i < size; ++i)
                ofs << "\t" << pDbl[i];
         } else {
            for (size_t i = 0; i < size; ++i)
                ofs << "\t" << pVar[i].value();
            ofs << "\nUncertainty";
            for (size_t i = 0; i < size; ++i)
                ofs << "\t" << pVar[i].uncertainty();
        }
        ofs << "\n" << EXPANSION_HEADER << "\n";
    }

    VarDbl value = outPrec? VarDbl(1, 0) : (isDbl? VarDbl(pDbl[0]) : pVar[0]);
    VarDbl variance = VarDbl();
    double var = this->variance();
    if (inPrec)
        var *= 1.0 /this->value() /this->value();
    double varn = var;
    VarDbl prevValue, prevVariance, newValue, newVariance;
    size_t monotonics = 0;
    size_t n = 2;
    for (; (n < _momentum.size()) && (n < size) && std::isfinite(varn) && (varn > 0); n += 2, varn *= var) {
        const VarDbl oldValue = value;
        const VarDbl oldVariance = variance;
        std::string except;
        size_t j = 1;
        try {
             if (isDbl)
                newValue = pDbl[n] * varn * _momentum[n];
            else
                newValue = pVar[n] * varn * _momentum[n];
            try {
                newVariance = VarDbl();
                if (isDbl) {
                    for (; j < n; ++j)
                        newVariance += pDbl[j] * pDbl[n - j] * varn * (_momentum[n] - _momentum[j] * _momentum[n - j]);
                } else {
                    for (; j < n; ++j)
                        newVariance += pVar[j] * pVar[n - j] * varn * (_momentum[n] - _momentum[j] * _momentum[n - j]);
                }
                try {
                    value += newValue;
                    try {
                        variance += newVariance;
                    } catch (InitException ex) {
                        except = "variance";
                    }
                } catch (InitException ex) {
                    except = "value";
                }
            } catch (InitException ex) {
                except = "newVariance";
            }
        } catch (InitException ex) {
            except = "newValue";
        }
        if (!except.empty()) {
            if (checkMonotonic && (monotonics >= MIN_MONOTONIC_COUNT)) { 
                value = oldValue;
                variance = oldVariance;
                break;
            }
            if (ofs) {
                ofs << "NotFiniteException\t" << except << "\t" << n << "\t" << varn << "\t" << _momentum[n] << "\t" << j << "\n";
                ofs.close();
            }
            throw NotFiniteException(name, *this, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }

        if (abs(newVariance.value()) <= abs(prevVariance.value()))
            monotonics += 1;
        else
            monotonics = 0;
        if (ofs) {
            ofs << n << "\t" << monotonics << "\t" << varn << "\t" << _momentum[n]
                << "\t" << value.value() << "\t" << value.uncertainty() 
                << "\t" << variance.value() << "\t" << variance.uncertainty()
                << "\t" << newValue.value() << "\t" << newValue.variance() 
                << "\t" << newVariance.value() << "\t" << newVariance.uncertainty() << "\n";
        }

        if (checkPositive && (variance.value() < 0)) {
            if (ofs) {
                ofs << "NotPositiveException\n";
                ofs.close();
            }
            throw NotPositiveException(name, *this, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }

        prevValue = newValue;
        prevVariance = newVariance;
    }

    if (checkMonotonic && (varn > 0) && (monotonics < MIN_MONOTONIC_COUNT)) {
        if (ofs) {
            ofs << "NotMonotonicException\n";
            ofs.close();
        }
        throw NotMonotonicException(name, *this, s1dTaylor, inPrec, outPrec,
                    value, variance, n, monotonics, newValue, newVariance);
    }
    const double unc = std::sqrt(variance.value()) *TAU;
    if (checkStability && !((std::abs(prevValue.value()) < unc) || (std::abs(prevValue.value()) < var_dbl::ulp(value.value())))) {
        if (ofs) {
            ofs << "NotStableException\t" << unc << "\t" << var_dbl::ulp(value.value()) << "\n";
            ofs.close();
        }
        throw NotStableException(name, *this, s1dTaylor, inPrec, outPrec,
                    value, variance, n, monotonics, newValue, newVariance);
    }

    if (outPrec) {
        try {
            if (isDbl) {
                value *= pDbl[0];
                variance *= pDbl[0] * pDbl[0];
            } else {
                value *= pVar[0];
                variance *= pVar[0] * pVar[0];
            }
        } catch (InitException ex) {
            if (ofs) {
                ofs << "NotFiniteException\tNormalization\n";
                ofs.close();
            }
            throw NotFiniteException(name, *this, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }
    }
    if (! std::isfinite(variance.value() + value.variance())) {
        if (ofs) {
            ofs << "NotFiniteException\tVariance overflow\n";
            ofs.close();
        }
        throw NotFiniteException(name, *this, s1dTaylor, inPrec, outPrec,
                    value, variance, n, monotonics, newValue, newVariance);
    }
    if (ofs) {
        ofs << OUTPUT_HEADER << "\n";
        ofs << value.value() << "\t" << value.variance() << "\t" << variance.value() << "\t" << variance.variance() << "\n";
        ofs.close();
    }
    return VarDbl(value.value(), variance.value() + value.variance(), true);
}

inline VarDbl VarDbl::polynominal(const std::vector<VarDbl>& sCoeff, const char* const dumpPath) const
{
    return polynominal(UnionVector(sCoeff), dumpPath);
}

inline VarDbl VarDbl::polynominal(const std::vector<double>& sCoeff, const char* const dumpPath) const
{
    return polynominal(UnionVector(sCoeff), dumpPath);
}


inline VarDbl VarDbl::polynominal(UnionVector sCoeff, const char* const dumpPath) const
{
    const bool isDbl = std::holds_alternative<std::vector<double>>(sCoeff);
    const size_t size = isDbl? std::get<std::vector<double>>(sCoeff).size() : std::get<std::vector<VarDbl>>(sCoeff).size();
    const double* const sDbl = isDbl? std::get<std::vector<double>>(sCoeff).data() : nullptr;
    const VarDbl* const sVar = isDbl? nullptr : std::get<std::vector<VarDbl>>(sCoeff).data();

    if (size >= VarDbl::_momentum.size()*2) {
        std::ostringstream oss;
        oss << "coefficient length " << size << " is too large!";
        throw std::invalid_argument(oss.str());
    }
    const int exp = size - 1;
    UnionVector s1dTaylor;
    if (isDbl)
        s1dTaylor = std::vector<double>(2*exp + 1);
    else
        s1dTaylor = std::vector<VarDbl>(2*exp + 1);
    double* const pDbl = isDbl? std::get<std::vector<double>>(s1dTaylor).data() : nullptr;
    VarDbl* const pVar = isDbl? nullptr : std::get<std::vector<VarDbl>>(s1dTaylor).data();
    if (isDbl)
        pDbl[0] = sDbl[0];
    else
        pVar[0] = sVar[0];
    std::vector<double> sPow;
    sPow.push_back(1);
    sPow.push_back(value());
    long sTaylor[2*exp + 1];
    sTaylor[0] = 1;
    for (int j = 1; j < size; ++j) {
        sTaylor[1] = j;
        for (int k = 2; k <= j; ++k) {
            sTaylor[k] = sTaylor[k - 1] * (j + 1 - k) / k;
            if (sPow.size() <= k)
                sPow.push_back(sPow.back() * value());
        }
        for (int k = 0; k <= j; ++k) {
            if (isDbl)
                pDbl[k] += sDbl[j] * sTaylor[k] * sPow[j - k];
            else
                pVar[k] += sVar[j] * sTaylor[k] * sPow[j - k];
        }
    }
    if (isDbl) {
        for (int j = 0; j < 2*exp + 1; ++j) {
            if (!std::isfinite(pDbl[j])) {
                std::ostringstream oss;
                oss << "polynominal coefficent at " << j << "/" << size;
                throw InitException(pDbl[j], 0, oss.str());
            }
        }
    } else {
        for (int j = 0; j < 2*exp + 1; ++j) {
            if (!std::isfinite(pVar[j].value()) || !std::isfinite(pVar[j].variance())) {
                std::ostringstream oss;
                oss << "polynominal coefficent at " << j << "/" << size;
                throw InitException(pVar[j].value(), pVar[j].variance(), oss.str());
            }
        }
    }
    std::ostringstream oss;
    oss << "Poly(" << exp << ")";
    return taylor1d(oss.str(), s1dTaylor, false, false, dumpPath,
                false, false, true, true);
}

inline VarDbl VarDbl::exp(const char* const dumpPath) const
{
    std::vector<double> sTaylor;
    sTaylor.reserve(_momentum.size());
    sTaylor.push_back(std::exp(value()));
    double n = 1;
    for (size_t i = 1; i < _momentum.size(); ++i) {
        n *= 1.0 / i;
        sTaylor.push_back(n);
    }
    std::ostringstream os;
    os << "exp(" << *this << ")";
    return taylor1d(os.str(), UnionVector(sTaylor), false, true, dumpPath);
}

inline VarDbl VarDbl::log(const char* const dumpPath) const 
{
    std::vector<double> sTaylor;
    sTaylor.reserve(_momentum.size());
    sTaylor.push_back(std::log(value()));
    for (size_t i = 1; i < _momentum.size(); ++i) {
        sTaylor.push_back((((i%2) == 1)? 1.0 : -1.0) / i);
    }
    std::ostringstream os;
    os << "log(" << *this << ")";
    return taylor1d(os.str(), UnionVector(sTaylor), true, false, dumpPath);
}

inline VarDbl VarDbl::sin(const char* const dumpPath) const 
{
    std::vector<double> sTaylor;
    sTaylor.reserve(_momentum.size());
    sTaylor.push_back(std::sin(value()));
    double n = 1;
    for (size_t i = 1; i < _momentum.size(); ++i) {
        n *= 1.0 / i;
        switch (i % 4) {
            case 0:
                sTaylor.push_back(std::sin(value()) * n);
                break;
            case 1:
                sTaylor.push_back(std::cos(value()) * n);  
                break;
            case 2:
                sTaylor.push_back(-std::sin(value()) * n);
                break;
            case 3:
                sTaylor.push_back(-std::cos(value()) * n);
                break;
        }
    }
    std::ostringstream os;
    os << "sin(" << *this << ")";
    return taylor1d(os.str(), UnionVector(sTaylor), false, false, dumpPath);
}

inline VarDbl VarDbl::pow(double exp, const char* const dumpPath) const 
{
    if (exp == 0)
        return VarDbl(1);
    else if (exp == 1)
        return *this;
    if ((exp > 0) && (std::ceil(exp) == std::floor(exp))) {
        std::vector<double> sCoeff(exp + 1);
        sCoeff[(size_t) exp] = 1;
        return polynominal(UnionVector(sCoeff));
    }
    std::vector<double> sTaylor;
    sTaylor.reserve(_momentum.size());
    sTaylor.push_back(std::pow(value(), exp));
    sTaylor.push_back(exp);
    for (size_t i = 2; i < _momentum.size(); ++i) {
        sTaylor.push_back(sTaylor[i - 1] * ((exp + 1)/i - 1));
    }
    std::ostringstream os;
    os << "pow(" << *this << ", " << exp << ")";
    return taylor1d(os.str(), UnionVector(sTaylor), true, true, dumpPath);
}

} // namespace var_dbl
#endif // __ValDbl_h__
