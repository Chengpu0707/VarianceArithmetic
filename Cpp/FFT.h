
#include <numbers>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "IndexSin.h"
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

    std::vector<VarDbl> transform(const std::vector<VarDbl>& sData, bool forward,
                                  bool traceSteps=false) const;
        // FFT of input {sData} of (real + imag) of size (2 << order)
        // When {traceSteps} is true, {ssStep} contains the data for intermediate steps,
        //  with [order + 1] for the result, and [order + 2] for the value error which could be non-exist
    template<typename T> requires std::floating_point<T> || std::integral<T>
    std::vector<VarDbl> transform(const std::vector<T>& sData, bool forward, 
                                  bool traceSteps=false) const;
    // intermediate steps for the transform() when {traceSteps} is true
    mutable std::vector<std::vector<VarDbl>> ssStep;

    FFT(IndexSin::SinSource sinSource = IndexSin::SinSource::Quart, const std::string& dumpDir = "") 
        : _sin(sinSource, dumpDir) {}
        // if {dumpDir} is not empty, read sin value from the file in the {dumpDir}

protected:
    const IndexSin _sin;
};



inline std::vector<size_t> FFT::bitReversedIndices(unsigned char order) 
{
    IndexSin::validateOrder(order);
    static std::unordered_map<unsigned char, std::vector<size_t>> _ssBitReversedIndex; 
    auto it = _ssBitReversedIndex.find(order);
    if (it != _ssBitReversedIndex.end())
        return it->second;
    std::vector<size_t> sRes(1 << order, 0);
    const size_t N = 1 << order;
    const size_t M = N >> 1;
    int i, j, k;
    j = 0;
    for (i = 0; i < N; i++) {
        sRes[i] = j;
        // next j from NumericalRecipesinC.pdf
        k = M; 
        while ((k != 0) && (j >= k)) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    _ssBitReversedIndex.insert({order, sRes});
    return sRes;
}

/*
    * 1-dimentional Fast Fourier Transformation (FFT)
    * 
    * @param sData     an array of size (2<<order), with each datum contains (real, image)
    * 
    * @return      an array of size (2<<order), with each datum contains (real, image)
    */
inline std::vector<VarDbl> FFT::transform(const std::vector<VarDbl>& sData, bool forward, bool traceSteps) const
{
    const unsigned char order = IndexSin::getOrder(sData.size() >> 1);
    const unsigned size = 1 << order;

    std::vector<VarDbl> sRes(2 << order);
    ssStep.clear();
    
    const std::vector<size_t> sIndex = bitReversedIndices(order);
    for (int i = 0; i < sIndex.size(); i++) {
        const int j = sIndex[i];
        sRes[(i << 1)] = sData[j << 1];
        sRes[(i << 1) + 1] = sData[(j << 1) + 1];
    }
    if (traceSteps)
        ssStep.push_back(sRes);

    for (int i = 0; i < (sIndex.size() - 1); i += 2 ) {
        const VarDbl rt = sRes[(i << 1)], it = sRes[(i << 1) + 1];
        sRes[(i << 1)] += sRes[(i << 1) + 2];
        sRes[(i << 1) + 1] += sRes[(i << 1) + 3];
        sRes[(i << 1) + 2] = rt - sRes[(i << 1) + 2];
        sRes[(i << 1) + 3] = it - sRes[(i << 1) + 3];
    }
    if (traceSteps)
        ssStep.push_back(sRes);

    for (unsigned o = 1, k = 4; o < order; ++o, k <<= 1) {
        for (long j = 0; j < (k >> 1); j++) {
            const VarDbl vcos = _sin.cos(j, o);
            const VarDbl vsin = _sin.sin(forward? j : -j, o);
            for (int i = 0; i < sIndex.size(); i += k ) {
                const int idx0 = (i + j) << 1;
                const int idx1 = idx0 + k;
                const VarDbl& r1 = sRes[idx1];
                const VarDbl& i1 = sRes[idx1 + 1];
        
                const VarDbl rd = r1 * vcos - i1 * vsin;
                const VarDbl id = i1 * vcos + r1 * vsin;
        
                sRes[idx1] = sRes[idx0] - rd;
                sRes[idx1 + 1] = sRes[idx0 + 1] - id;
                sRes[idx0] += rd;
                sRes[idx0 + 1] += id;
            }   // for( i
        }
        if (traceSteps)
            ssStep.push_back(sRes);
    }
    
    if (!forward) {
        for (int i = 0; i < (sIndex.size() << 1); i ++ ) {
            sRes[i] *= 1.0/size;
        }
    }
    if (traceSteps)
        ssStep.push_back(sRes);
    return sRes;
}

template<typename T> requires std::floating_point<T> || std::integral<T>
inline std::vector<VarDbl> FFT::transform(const std::vector<T>& sData, bool forward, bool traceSteps) const
{
    return transform(std::vector<VarDbl>{sData}, forward, traceSteps);
}


} // namespace var_dbl
#endif  // __FFT_h__