#include <exception>
#if __cplusplus >= 201103L
#include <functional>
#endif
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#if __cplusplus >= 201703L
#include <variant>
#endif
#include <vector>

#include "Momentum.h"
#include "VarDbl.h"

#ifndef __Taylor_h__
#define __Taylor_h__
namespace var_dbl
{

#if __cplusplus >= 201703L
using UnionVector = std::variant<std::vector<VarDbl>, std::vector<double>>;
#else
struct UnionVector {
    std::vector<VarDbl> var;
    std::vector<double> dbl;
    bool isDbl;
    UnionVector() : isDbl(true) {}
    UnionVector(const std::vector<double>& v) : dbl(v), isDbl(true) {}
    UnionVector(const std::vector<VarDbl>& v) : var(v), isDbl(false) {}
#if __cplusplus >= 201103L
    UnionVector(std::vector<double>&& v) : dbl(std::move(v)), isDbl(true) {}
    UnionVector(std::vector<VarDbl>&& v) : var(std::move(v)), isDbl(false) {}
#endif
};
#endif

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

#if __cplusplus < 201103L
    virtual ~TaylorIdException() throw() {}
#endif
    virtual std::string name() const { return "TaylorIdException"; }
};

struct NotFiniteException : public TaylorIdException
{
#if __cplusplus >= 201103L
    using TaylorIdException::TaylorIdException;
#else
    explicit NotFiniteException(const std::string what,
                const VarDbl& input, const UnionVector& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance,
                size_t n, size_t monotonics, const VarDbl& newValue, const VarDbl& newVariance) :
            TaylorIdException(what, input, s1dTaylor, inPrec, outPrec,
                value, variance, n, monotonics, newValue, newVariance) {}
    virtual ~NotFiniteException() throw() {}
#endif
    virtual std::string name() const { return "NotFiniteException"; }
};

struct NotReliableException : public TaylorIdException
{
#if __cplusplus >= 201103L
    using TaylorIdException::TaylorIdException;
#else
    explicit NotReliableException(const std::string what,
                const VarDbl& input, const UnionVector& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance,
                size_t n, size_t monotonics, const VarDbl& newValue, const VarDbl& newVariance) :
            TaylorIdException(what, input, s1dTaylor, inPrec, outPrec,
                value, variance, n, monotonics, newValue, newVariance) {}
    virtual ~NotReliableException() throw() {}
#endif
    virtual std::string name() const { return "NotReliableException"; }
};

struct NotMonotonicException : public TaylorIdException
{
#if __cplusplus >= 201103L
    using TaylorIdException::TaylorIdException;
#else
    explicit NotMonotonicException(const std::string what,
                const VarDbl& input, const UnionVector& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance,
                size_t n, size_t monotonics, const VarDbl& newValue, const VarDbl& newVariance) :
            TaylorIdException(what, input, s1dTaylor, inPrec, outPrec,
                value, variance, n, monotonics, newValue, newVariance) {}
    virtual ~NotMonotonicException() throw() {}
#endif
    virtual std::string name() const { return "NotMonotonicException"; }
};

struct NotStableException : public TaylorIdException
{
#if __cplusplus >= 201103L
    using TaylorIdException::TaylorIdException;
#else
    explicit NotStableException(const std::string what,
                const VarDbl& input, const UnionVector& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance,
                size_t n, size_t monotonics, const VarDbl& newValue, const VarDbl& newVariance) :
            TaylorIdException(what, input, s1dTaylor, inPrec, outPrec,
                value, variance, n, monotonics, newValue, newVariance) {}
    virtual ~NotStableException() throw() {}
#endif
    virtual std::string name() const { return "NotStableException"; }
};

struct NotPositiveException : public TaylorIdException
{
#if __cplusplus >= 201103L
    using TaylorIdException::TaylorIdException;
#else
    explicit NotPositiveException(const std::string what,
                const VarDbl& input, const UnionVector& s1dTaylor, bool inPrec, bool outPrec,
                const VarDbl& value, const VarDbl& variance,
                size_t n, size_t monotonics, const VarDbl& newValue, const VarDbl& newVariance) :
            TaylorIdException(what, input, s1dTaylor, inPrec, outPrec,
                value, variance, n, monotonics, newValue, newVariance) {}
    virtual ~NotPositiveException() throw() {}
#endif
    virtual std::string name() const { return "NotPositiveException"; }
};


const NormalMomentum IDEAL_MOMENTUM(5);


/*
A static class for Taylor expansion.
It is not intended to be instantiated, so the constructor is deleted.
*/
class Taylor {
#if __cplusplus >= 201103L
    Taylor() = delete;
    Taylor(const Taylor&) = delete;
    Taylor(Taylor&&) = delete;
    Taylor& operator=(Taylor&) = delete;
    Taylor& operator=(Taylor&&) = delete;
    Taylor& operator=(const Taylor&) = delete;
#else
private:
    Taylor();
    Taylor(const Taylor&);
    Taylor& operator=(const Taylor&);
public:
#endif

#if __cplusplus >= 201103L
    constexpr static const int MIN_MONOTONIC_COUNT = 20;
    constexpr static const char* INPUT_HEADER =
            "result\tvalue\tuncertainty\tinPrec\toutPrec\tbounding\tmaxOrder\tMinMonotonic\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive\tName";
    constexpr static const char* EXPANSION_HEADER =
            "Order\tTaylor Value\tTaylor Uncertainty\tExponent\tMomentum\tMonotonics\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty";
    constexpr static const char* OUTPUT_HEADER =
            "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tException";
#else
    static const int MIN_MONOTONIC_COUNT;
    static const char* INPUT_HEADER;
    static const char* EXPANSION_HEADER;
    static const char* OUTPUT_HEADER;
#endif

    static VarDbl taylor1d(const VarDbl& in, std::string name, const UnionVector s1dTaylor, bool inPrec, bool outPrec,
            const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM,
            bool checkMonotonic=true, bool checkStability=true, bool checkReliablity=true, bool checkPositive=true);

    static VarDbl poly1d(const VarDbl& in, UnionVector sCoeff,
            const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);

public:
    static VarDbl exp(const VarDbl& in,
            const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
    static VarDbl log(const VarDbl& in,
            const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
    static VarDbl sin(const VarDbl& in,
            const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
    static VarDbl pow(const VarDbl& in, double exp,
        const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
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
            momentum.leakage-fold of the expansion uncertainty.
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
    static VarDbl taylor1d(const VarDbl& in, std::string name, const std::vector<double>& s1dTaylor, bool inPrec, bool outPrec,
#if __cplusplus >= 201103L
                    const char* const dumpPath=nullptr, const Momentum& momentum=IDEAL_MOMENTUM,
#else
                    const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM,
#endif
                    bool checkMonotonic=true, bool checkStability=true, bool checkReliablity=true, bool checkPositive=true);
    static VarDbl taylor1d(const VarDbl& in, std::string name, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
#if __cplusplus >= 201103L
                    const char* const dumpPath=nullptr, const Momentum& momentum=IDEAL_MOMENTUM,
#else
                    const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM,
#endif
                    bool checkMonotonic=true, bool checkStability=true, bool checkReliablity=true, bool checkPositive=true);

    /*
     1d polynomial with coeeficient {sCoeff}
    */
                    static VarDbl poly1d(const VarDbl& in, const std::vector<double>& sCoeff,
                       const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
    static VarDbl poly1d(const VarDbl& in, const std::vector<VarDbl>& sCoeff,
                       const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
#if __cplusplus >= 202002L
    template <typename T> requires std::integral<T>
#else
    template <typename T>
#endif
    static VarDbl poly1d(const VarDbl& in, const std::vector<T>& sCoeff,
                       const char* const dumpPath=0, const Momentum& momentum=IDEAL_MOMENTUM);
};


#if __cplusplus < 201103L
const int Taylor::MIN_MONOTONIC_COUNT = 20;
const char* Taylor::INPUT_HEADER =
        "result\tvalue\tuncertainty\tinPrec\toutPrec\tbounding\tmaxOrder\tMinMonotonic\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive\tName";
const char* Taylor::EXPANSION_HEADER =
        "Order\tTaylor Value\tTaylor Uncertainty\tExponent\tMomentum\tMonotonics\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty";
const char* Taylor::OUTPUT_HEADER =
        "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tException";
#endif


inline VarDbl Taylor::taylor1d(
    const VarDbl& in, std::string name, const std::vector<double>& s1dTaylor, bool inPrec, bool outPrec,
    const char* const dumpPath, const Momentum& momentum,
    bool checkMonotonic, bool checkStability, bool checkReliablity, bool checkPositive
    )
{
#if __cplusplus >= 201103L
    return taylor1d(in, name, UnionVector(std::move(s1dTaylor)), inPrec, outPrec,
                    dumpPath, momentum,
                    checkMonotonic, checkStability, checkReliablity, checkPositive);
#else
    return taylor1d(in, name, UnionVector(s1dTaylor), inPrec, outPrec,
                    dumpPath, momentum,
                    checkMonotonic, checkStability, checkReliablity, checkPositive);
#endif
}


inline VarDbl Taylor::taylor1d(
    const VarDbl& in, std::string name, const std::vector<VarDbl>& s1dTaylor, bool inPrec, bool outPrec,
    const char* const dumpPath, const Momentum& momentum,
    bool checkMonotonic, bool checkStability, bool checkReliablity, bool checkPositive
    )
{
#if __cplusplus >= 201103L
    return taylor1d(in, name, UnionVector(std::move(s1dTaylor)), inPrec, outPrec,
                    dumpPath, momentum,
                    checkMonotonic, checkStability, checkReliablity, checkPositive);
#else
    return taylor1d(in, name, UnionVector(s1dTaylor), inPrec, outPrec,
                    dumpPath, momentum,
                    checkMonotonic, checkStability, checkReliablity, checkPositive);
#endif
}


inline VarDbl Taylor::taylor1d(
    const VarDbl& in, std::string name, const UnionVector s1dTaylor, bool inPrec, bool outPrec,
    const char* const dumpPath, const Momentum& momentum,
    bool checkMonotonic, bool checkStability, bool checkReliablity, bool checkPositive
    )
{
#if __cplusplus >= 201703L
    const bool isDbl = std::holds_alternative<std::vector<double>>(s1dTaylor);
    const size_t size = isDbl? std::get<std::vector<double>>(s1dTaylor).size() : std::get<std::vector<VarDbl>>(s1dTaylor).size();
    const double* pDbl = isDbl? std::get<std::vector<double>>(s1dTaylor).data() : nullptr;
    const VarDbl* pVar = isDbl? nullptr : std::get<std::vector<VarDbl>>(s1dTaylor).data();
#else
    const bool isDbl = s1dTaylor.isDbl;
    const size_t size = isDbl? s1dTaylor.dbl.size() : s1dTaylor.var.size();
    const double* pDbl = isDbl? s1dTaylor.dbl.data() : NULL;
    const VarDbl* pVar = isDbl? NULL : s1dTaylor.var.data();
#endif

    if (isDbl) {
        for (int j = 0; j < (int)size; ++j) {
            if (!std::isfinite(pDbl[j])) {
                std::ostringstream oss;
                oss << "Taylor " << pDbl[j] << " at " << j << "/" << size;
                throw std::invalid_argument(oss.str());
            }
        }
    } else {
        for (int j = 0; j < (int)size; ++j) {
            if (!std::isfinite(pVar[j].value()) || !std::isfinite(pVar[j].variance())) {
                std::ostringstream oss;
                oss << "Taylor " << pVar[j] << " at " << j << "/" << size;
                throw std::invalid_argument(oss.str());
            }
        }
    }
    if (in.uncertainty() == 0) {
        if (isDbl)
            return VarDbl(pDbl[0]);
        else
            return pVar[0];
    }
    std::ofstream ofs(dumpPath);
    if (ofs) {
        ofs << std::setprecision(15);
        ofs << INPUT_HEADER <<"\n";
        ofs << (isDbl? pDbl[0] : pVar[0].value()) << '\t' << in.value() << '\t' << in.uncertainty() << '\t'
            << inPrec << '\t' << outPrec << '\t'
            << momentum.bounding << '\t' << momentum.maxOrder() << '\t' <<Taylor::MIN_MONOTONIC_COUNT << '\t'
            << checkMonotonic << '\t' << checkStability << '\t' << checkReliablity << '\t' << checkPositive << '\t'
            << name << '\n';
        ofs << EXPANSION_HEADER << "\n";
    }

    VarDbl value = outPrec? VarDbl(1, 0) : (isDbl? VarDbl(pDbl[0]) : pVar[0]);
    VarDbl variance, newValue, newVariance, prevVariance;
#if __cplusplus >= 201103L
    auto dumpOutput = [&ofs](const VarDbl& value, const VarDbl& variance, const std::string exception) {
        if (ofs) {
            ofs << OUTPUT_HEADER <<"\n";
            ofs << value.value() << '\t' << value.uncertainty() << '\t' << variance.value() << '\t' << value.uncertainty() << '\t'
                << exception << "\n";
            ofs.close();
        }
    };
#else
    struct _DumpOutput {
        std::ofstream& ofs;
        _DumpOutput(std::ofstream& ofs_) : ofs(ofs_) {}
        void operator()(const VarDbl& value, const VarDbl& variance, const std::string& exception) const {
            if (ofs) {
                ofs << OUTPUT_HEADER <<"\n";
                ofs << value.value() << '\t' << value.uncertainty() << '\t' << variance.value() << '\t' << value.uncertainty() << '\t'
                    << exception << "\n";
                const_cast<std::ofstream&>(ofs).close();
            }
        }
    } dumpOutput(ofs);
#endif

    const double unc = inPrec? in.uncertainty() /in.value() : in.uncertainty();
    double uncN = unc;
    size_t monotonics = 0;
    bool monotonicPrev = true, finite=true;
    size_t n = 1;
    for (; (n < momentum.maxOrder()) && (n < size) && std::isfinite(uncN) && (uncN > 0); ++n, uncN *= unc) {
        std::string except;
        size_t j = 1;
        try {
            if (isDbl) {
                if (abs(uncN) < 1)
                    newValue = pDbl[n] * (uncN * momentum[n]);
                else
                    newValue = pDbl[n] * uncN * momentum[n];
            }
            else {
                if (abs(uncN) < 1)
                    newValue = pVar[n] * (uncN * momentum[n]);
                else
                    newValue = pVar[n] * uncN * momentum[n];
            }
            try {
                newVariance = VarDbl();
                if (isDbl) {
                    double newVar = 0;
                    for (; j < n; ++j) {
                        if ((uncN < 1))
                            newVar += pDbl[j] * pDbl[n - j] * (uncN * (momentum[n] - momentum[j] * momentum[n - j]));
                        else
                            newVar += pDbl[j] * pDbl[n - j] * uncN * (momentum[n] - momentum[j] * momentum[n - j]);
                    }
                    newVariance += newVar;
                } else {
                    for (; j < n; ++j) {
                        if ((uncN < 1))
                            newVariance += pVar[j] * pVar[n - j] * (uncN * (momentum[n] - momentum[j] * momentum[n - j]));
                        else
                            newVariance += pVar[j] * pVar[n - j] * uncN * (momentum[n] - momentum[j] * momentum[n - j]);
                    }
                }
                try {
                    value += newValue;
                    try {
                        variance += newVariance;
                    } catch (InitException ex) {
                        except = "NotFiniteException\tvariance";
                    }
                } catch (InitException ex) {
                    except = "NotFiniteException\tvalue";
                }
            } catch (InitException ex) {
                except = "NotFiniteException\tnewVariance";
            }
        } catch (InitException ex) {
            except = "NotFiniteException\tnewValue";
        }


        if (except.empty())
        {
            if (!std::isfinite(newValue.value()) || !std::isfinite(newValue.uncertainty()))
                except = "NotFiniteException newValue";
            else if (!std::isfinite(newVariance.value()) || !std::isfinite(newVariance.uncertainty()))
                except = "NotFiniteException newValue";
            else if (!std::isfinite(value.value()) || !std::isfinite(value.uncertainty()))
                except = "NotFiniteException value";
            else if (!std::isfinite(variance.value()) || !std::isfinite(variance.uncertainty()))
                except = "NotFiniteException newValue";
        }
        if (!except.empty()) {
            dumpOutput(value, variance, except);
            throw NotFiniteException(name, in, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }

        if ((n & 1) == 0) {
            if (abs(newVariance.value()) <= abs(prevVariance.value()))
                monotonics += 1;
            else if ((monotonics >= Taylor::MIN_MONOTONIC_COUNT) && monotonicPrev)
                monotonicPrev = false;
            else
                monotonics = 0;
            prevVariance = newVariance;
        }
        if (ofs) {
            ofs << n
                << "\t" << (isDbl? pDbl[n] : pVar[n].value()) << "\t" << (isDbl? 0 : pVar[n].uncertainty())
                << "\t" << uncN << "\t" << momentum[n] << "\t" << monotonics
                << "\t" << value.value() << "\t" << value.uncertainty()
                << "\t" << variance.value() << "\t" << variance.uncertainty()
                << "\t" << newValue.value() << "\t" << newValue.variance()
                << "\t" << newVariance.value() << "\t" << newVariance.uncertainty() << "\n";
        }

        if (checkPositive && (variance.value() < 0)) {
            dumpOutput(value, variance, "NotPositiveException");
            throw NotPositiveException(name, in, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }
    }

    if (checkMonotonic && (uncN > 0) && (monotonics < MIN_MONOTONIC_COUNT)) {
        dumpOutput(value, variance, "NotMonotonicException");
        throw NotMonotonicException(name, in, s1dTaylor, inPrec, outPrec,
                    value, variance, n, monotonics, newValue, newVariance);
    }
    if (checkStability) {
        if (std::abs(newValue.value()) > var_dbl::ulp(value.value())) {
            dumpOutput(value, variance, "NotStableException\tValue");
            throw NotStableException(name, in, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }
        if (std::abs(newValue.value()) > std::sqrt(variance.value() + value.variance()) *momentum.leakage) {
            dumpOutput(value, variance, "NotStableException\tVariance");
            throw NotStableException(name, in, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }
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
            dumpOutput(value, variance, "NotFiniteException\tNormalization");
            throw NotFiniteException(name, in, s1dTaylor, inPrec, outPrec,
                        value, variance, n, monotonics, newValue, newVariance);
        }
    }
    if (! std::isfinite(variance.value() + value.variance())) {
        dumpOutput(value, variance, "NotFiniteException\tVariance Overflow");
        throw NotFiniteException(name, in, s1dTaylor, inPrec, outPrec,
                    value, variance, n, monotonics, newValue, newVariance);
    }
    dumpOutput(value, variance, "");
    return VarDbl(value.value(), std::sqrt(variance.value() + value.variance()));
}

inline VarDbl Taylor::poly1d(const VarDbl& in, const std::vector<VarDbl>& sCoeff,
        const char* const dumpPath, const Momentum& momentum)
{
    return poly1d(in, UnionVector(sCoeff), dumpPath, momentum);
}

inline VarDbl Taylor::poly1d(const VarDbl& in, const std::vector<double>& sCoeff,
        const char* const dumpPath, const Momentum& momentum)
{
    return poly1d(in, UnionVector(sCoeff), dumpPath, momentum);
}

#if __cplusplus >= 202002L
template <typename T> requires std::integral<T>
#else
template <typename T>
#endif
inline VarDbl Taylor::poly1d(const VarDbl& in, const std::vector<T>& sCoeff,
        const char* const dumpPath, const Momentum& momentum)
{
    std::vector<double> inCoeff(sCoeff.size());
    for (size_t i = 0; i < sCoeff.size(); ++i) {
        inCoeff[i] = static_cast<double>(sCoeff[i]);
    }
    return poly1d(in, UnionVector(inCoeff), dumpPath, momentum);
}


inline VarDbl Taylor::poly1d(const VarDbl& in, UnionVector sCoeff,
        const char* const dumpPath, const Momentum& momentum)
{
#if __cplusplus >= 201703L
    const bool isDbl = std::holds_alternative<std::vector<double>>(sCoeff);
    const size_t size = isDbl? std::get<std::vector<double>>(sCoeff).size() : std::get<std::vector<VarDbl>>(sCoeff).size();
    const double* const sDbl = isDbl? std::get<std::vector<double>>(sCoeff).data() : nullptr;
    const VarDbl* const sVar = isDbl? nullptr : std::get<std::vector<VarDbl>>(sCoeff).data();
#else
    const bool isDbl = sCoeff.isDbl;
    const size_t size = isDbl? sCoeff.dbl.size() : sCoeff.var.size();
    const double* const sDbl = isDbl? sCoeff.dbl.data() : NULL;
    const VarDbl* const sVar = isDbl? NULL : sCoeff.var.data();
#endif

    if (size >= momentum.maxOrder()*2) {
        std::ostringstream oss;
        oss << "coefficient length " << size << " is too large!";
        throw std::invalid_argument(oss.str());
    }
    const int exp = (int)(size - 1);
    UnionVector s1dTaylor;
#if __cplusplus >= 201703L
    if (isDbl)
        s1dTaylor = std::vector<double>(2*exp + 1);
    else
        s1dTaylor = std::vector<VarDbl>(2*exp + 1);
    double* const pDbl = isDbl? std::get<std::vector<double>>(s1dTaylor).data() : nullptr;
    VarDbl* const pVar = isDbl? nullptr : std::get<std::vector<VarDbl>>(s1dTaylor).data();
#else
    if (isDbl) {
        s1dTaylor.dbl.assign(2*exp + 1, 0.0);
        s1dTaylor.var.clear();
        s1dTaylor.isDbl = true;
    } else {
        s1dTaylor.var.assign(2*exp + 1, VarDbl());
        s1dTaylor.dbl.clear();
        s1dTaylor.isDbl = false;
    }
    double* const pDbl = isDbl? s1dTaylor.dbl.data() : NULL;
    VarDbl* const pVar = isDbl? NULL : s1dTaylor.var.data();
#endif
    if (isDbl)
        pDbl[0] = sDbl[0];
    else
        pVar[0] = sVar[0];
    std::vector<double> sPow;
    sPow.push_back(1);
    sPow.push_back(in.value());
    long sTaylor[2*exp + 1];
    sTaylor[0] = 1;
    for (int j = 1; j < (int)size; ++j) {
        sTaylor[1] = j;
        for (int k = 2; k <= j; ++k) {
            sTaylor[k] = sTaylor[k - 1] * (j + 1 - k) / k;
            if ((int)sPow.size() <= k)
                sPow.push_back(sPow.back() * in.value());
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
                oss << "poly1d coefficent at " << j << "/" << size;
                throw InitException(pDbl[j], 0, oss.str());
            }
        }
    } else {
        for (int j = 0; j < 2*exp + 1; ++j) {
            if (!std::isfinite(pVar[j].value()) || !std::isfinite(pVar[j].variance())) {
                std::ostringstream oss;
                oss << "poly1d coefficent at " << j << "/" << size;
                throw InitException(pVar[j].value(), pVar[j].variance(), oss.str());
            }
        }
    }
    std::ostringstream oss;
    oss << "Poly(" << exp << ")";
    return taylor1d(in, oss.str(), s1dTaylor, false, false, dumpPath, momentum,
                false, false, true, true);
}

inline VarDbl Taylor::exp(const VarDbl& in,
        const char* const dumpPath, const Momentum& momentum)
{
    std::vector<double> sTaylor;
    sTaylor.reserve(momentum.maxOrder());
    sTaylor.push_back(std::exp(in.value()));
    double n = 1;
    for (size_t i = 1; i < momentum.maxOrder(); ++i) {
        n *= 1.0 / i;
        sTaylor.push_back(n);
    }
    std::ostringstream os;
    os << "exp(" << in << ")";
    return taylor1d(in, os.str(), UnionVector(sTaylor), false, true, dumpPath, momentum);
}

inline VarDbl Taylor::log(const VarDbl& in,
        const char* const dumpPath, const Momentum& momentum)
{
    std::vector<double> sTaylor;
    sTaylor.reserve(momentum.maxOrder());
    sTaylor.push_back(std::log(in.value()));
    for (size_t i = 1; i < momentum.maxOrder(); ++i) {
        sTaylor.push_back((((i%2) == 1)? 1.0 : -1.0) / i);
    }
    std::ostringstream os;
    os << "log(" << in << ")";
    return taylor1d(in, os.str(), UnionVector(sTaylor), true, false, dumpPath, momentum);
}

inline VarDbl Taylor::sin(const VarDbl& in,
        const char* const dumpPath, const Momentum& momentum)
{
    std::vector<double> sTaylor;
    sTaylor.reserve(momentum.maxOrder());
    sTaylor.push_back(std::sin(in.value()));
    double n = 1;
    for (size_t i = 1; i < momentum.maxOrder(); ++i) {
        n *= 1.0 / i;
        switch (i % 4) {
            case 0:
                sTaylor.push_back(std::sin(in.value()) * n);
                break;
            case 1:
                sTaylor.push_back(std::cos(in.value()) * n);
                break;
            case 2:
                sTaylor.push_back(-std::sin(in.value()) * n);
                break;
            case 3:
                sTaylor.push_back(-std::cos(in.value()) * n);
                break;
        }
    }
    std::ostringstream os;
    os << "sin(" << in << ")";
    return taylor1d(in, os.str(), UnionVector(sTaylor), false, false, dumpPath, momentum);
}

inline VarDbl Taylor::pow(const VarDbl& in, double exp,
        const char* const dumpPath, const Momentum& momentum)
{
    if (exp == 0)
        return VarDbl(1);
    else if (exp == 1)
        return in;
    bool neg = in.value() < 0;
    if (std::ceil(exp) == std::floor(exp)) {
        if (exp > 0) {
            std::vector<double> sCoeff((size_t)(exp + 1));
            sCoeff[(size_t) exp] = 1;
            return poly1d(in, UnionVector(sCoeff), dumpPath, momentum);
        }
        if ((((int) exp) & 1) == 0)
            neg = false;
    } else if (neg) {
        throw std::invalid_argument("Negative base with non-integer exponent");
    }
    std::vector<VarDbl> sTaylor;
    sTaylor.reserve(momentum.maxOrder());
    sTaylor.push_back(VarDbl(std::pow(std::abs(in.value()), exp)));
    sTaylor.push_back(VarDbl(exp));
    for (size_t i = 2; i < momentum.maxOrder(); ++i) {
        sTaylor.push_back(sTaylor[i - 1] * ((exp + 1)/i - 1));
    }
    std::ostringstream os;
    os << "pow(" << in << ", " << exp << ")";
    VarDbl res = taylor1d((in.value() >= 0)? in : -in,
                          os.str(), UnionVector(sTaylor), true, true, dumpPath, momentum);
    return neg? -res : res;
}

inline VarDbl VarDbl::operator/=(const VarDbl& other)
{
    return *this *=  Taylor::pow(other, -1);
}

} // namespace var_dbl
#endif // __Taylor_h__
