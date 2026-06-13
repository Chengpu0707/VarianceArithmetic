/*
Variance-arithmetic FFT: forward/inverse Fast Fourier Transform on VarDbl
sequences using IndexSin for trig accuracy. Tracks bit-reversed indices and
operates with statistical Taylor expansion to propagate uncertainty through
each butterfly stage.
*/
#if __cplusplus >= 201103L
#include <unordered_map>
#else
#include <map>
#endif
#if __cplusplus >= 201703L
#include <type_traits>
#endif
#include <sstream>
#include <vector>

#include "IndexSin.h"
#include "Interval.h"
#include "Taylor.h"
#include "Stat.h"

#ifndef __FFT_h__
#define __FFT_h__
namespace var_dbl
{

/*
 * A class to facilitate testing of FFT
 */
class FFT {
public:
    static std::vector<size_t> bitReversedIndices(unsigned char order);

    // Templated FFT over T in {VarDbl, Interval}.  Same butterfly body for both;
    // traceSteps is honored only when T == VarDbl (ssStep is typed to VarDbl).
    // Input/output are size (2<<order) arrays of interleaved (real, imag) pairs.
#if __cplusplus >= 202002L
    template<typename T>
        requires (std::is_same_v<T, VarDbl> || std::is_same_v<T, Interval>)
    std::vector<T> transform(const std::vector<T>& sData, bool forward,
                              bool traceSteps=false) const;
#elif __cplusplus >= 201103L
    template<typename T,
             typename = typename std::enable_if<
                std::is_same<T, VarDbl>::value || std::is_same<T, Interval>::value>::type>
    std::vector<T> transform(const std::vector<T>& sData, bool forward,
                              bool traceSteps=false) const;
#endif

    // Promote arithmetic-typed inputs (int, double, ...) to VarDbl first.
#if __cplusplus >= 202002L
    template<typename T> requires std::floating_point<T> || std::integral<T>
    std::vector<VarDbl> transform(const std::vector<T>& sData, bool forward,
                                  bool traceSteps=false) const;
#elif __cplusplus >= 201103L
    template<typename T,
             typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    std::vector<VarDbl> transform(const std::vector<T>& sData, bool forward,
                                  bool traceSteps=false) const;
#endif

    // intermediate steps for the transform() when {traceSteps} is true (VarDbl only)
    mutable std::vector<std::vector<VarDbl> > ssStep;

    FFT(IndexSin::SinSource sinSource = IndexSin::Quart, const std::string& dumpDir = "")
        : _sin(sinSource, dumpDir) {}
        // if {dumpDir} is not empty, read sin value from the file in the {dumpDir}

protected:
    const IndexSin _sin;
};



inline std::vector<size_t> FFT::bitReversedIndices(unsigned char order)
{
    IndexSin::validateOrder(order);
#if __cplusplus >= 201103L
    static std::unordered_map<unsigned char, std::vector<size_t>> _ssBitReversedIndex;
    auto it = _ssBitReversedIndex.find(order);
#else
    static std::map<unsigned char, std::vector<size_t> > _ssBitReversedIndex;
    std::map<unsigned char, std::vector<size_t> >::iterator it = _ssBitReversedIndex.find(order);
#endif
    if (it != _ssBitReversedIndex.end())
        return it->second;
    std::vector<size_t> sRes(1 << order, 0);
    const size_t N = 1 << order;
    const size_t M = N >> 1;
    int i, j, k;
    j = 0;
    for (i = 0; i < (int)N; i++) {
        sRes[i] = j;
        // next j from NumericalRecipesinC.pdf
        k = M;
        while ((k != 0) && (j >= k)) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
#if __cplusplus >= 201103L
    _ssBitReversedIndex.insert({order, sRes});
#else
    _ssBitReversedIndex.insert(std::make_pair(order, sRes));
#endif
    return sRes;
}

/*
 * 1-dimensional Fast Fourier Transformation (FFT) over T in {VarDbl, Interval}.
 *
 * @param sData     size (2<<order) array of interleaved (real, imag) pairs.
 * @param forward   true: forward DFT; false: inverse DFT (with 1/N normalization).
 * @param traceSteps  honored only when T == VarDbl; populates ssStep with the
 *                  intermediate state after each butterfly stage.
 *
 * Twiddles come from IndexSin::cos/sin (VarDbl-valued) and are converted via the
 * one-argument constructor T(VarDbl) — copy for VarDbl, [v.value-w, v.value+w]
 * (with w = max(unc, ulp(value))) for Interval.
 */
#if __cplusplus >= 202002L
template<typename T>
    requires (std::is_same_v<T, VarDbl> || std::is_same_v<T, Interval>)
inline std::vector<T> FFT::transform(const std::vector<T>& sData, bool forward, bool traceSteps) const
#elif __cplusplus >= 201103L
template<typename T, typename>
inline std::vector<T> FFT::transform(const std::vector<T>& sData, bool forward, bool traceSteps) const
#endif
{
    const unsigned char order = IndexSin::getOrder(sData.size() >> 1);
    const unsigned size = 1 << order;

    std::vector<T> sRes(2 << order);
#if __cplusplus >= 201703L
    if constexpr (std::is_same_v<T, VarDbl>) {
        ssStep.clear();
        if (traceSteps) ssStep.push_back(sData);
    }
#endif

    const std::vector<size_t> sIndex = bitReversedIndices(order);
    for (size_t i = 0; i < sIndex.size(); ++i) {
        const size_t j = sIndex[i];
        sRes[(i << 1)]     = sData[j << 1];
        sRes[(i << 1) + 1] = sData[(j << 1) + 1];
    }
#if __cplusplus >= 201703L
    if constexpr (std::is_same_v<T, VarDbl>) {
        if (traceSteps) ssStep.push_back(sRes);
    }
#endif

    for (size_t i = 0; i + 1 < sIndex.size(); i += 2) {
        const T rt = sRes[(i << 1)];
        const T it = sRes[(i << 1) + 1];
        sRes[(i << 1)]     += sRes[(i << 1) + 2];
        sRes[(i << 1) + 1] += sRes[(i << 1) + 3];
        sRes[(i << 1) + 2] = rt - sRes[(i << 1) + 2];
        sRes[(i << 1) + 3] = it - sRes[(i << 1) + 3];
    }
#if __cplusplus >= 201703L
    if constexpr (std::is_same_v<T, VarDbl>) {
        if (traceSteps) ssStep.push_back(sRes);
    }
#endif

    for (unsigned o = 1, k = 4; o < order; ++o, k <<= 1) {
        for (long j = 0; j < (long)(k >> 1); ++j) {
            const T vcos(_sin.cos(j, o));
            const T vsin(_sin.sin(forward ? j : -j, o));
            for (size_t i = 0; i < sIndex.size(); i += k) {
                const size_t idx0 = (i + j) << 1;
                const size_t idx1 = idx0 + k;
                const T r1 = sRes[idx1];
                const T i1 = sRes[idx1 + 1];
                const T rd = r1 * vcos - i1 * vsin;
                const T id = i1 * vcos + r1 * vsin;
                sRes[idx1]     = sRes[idx0]     - rd;
                sRes[idx1 + 1] = sRes[idx0 + 1] - id;
                sRes[idx0]     += rd;
                sRes[idx0 + 1] += id;
            }
        }
#if __cplusplus >= 201703L
        if constexpr (std::is_same_v<T, VarDbl>) {
            if (traceSteps) ssStep.push_back(sRes);
        }
#endif
    }

    if (!forward) {
        const double invN = 1.0 / size;
        for (size_t i = 0; i < (sIndex.size() << 1); ++i)
            sRes[i] *= invN;
    }
#if __cplusplus >= 201703L
    if constexpr (std::is_same_v<T, VarDbl>) {
        if (traceSteps) ssStep.push_back(sRes);
    }
#endif
    return sRes;
}

#if __cplusplus >= 202002L
template<typename T> requires std::floating_point<T> || std::integral<T>
inline std::vector<VarDbl> FFT::transform(const std::vector<T>& sData, bool forward, bool traceSteps) const
{
    return transform<VarDbl>(std::vector<VarDbl>(sData.begin(), sData.end()), forward, traceSteps);
}
#elif __cplusplus >= 201103L
template<typename T, typename>
inline std::vector<VarDbl> FFT::transform(const std::vector<T>& sData, bool forward, bool traceSteps) const
{
    return transform<VarDbl>(std::vector<VarDbl>(sData.begin(), sData.end()), forward, traceSteps);
}
#endif


} // namespace var_dbl
#endif  // __FFT_h__
