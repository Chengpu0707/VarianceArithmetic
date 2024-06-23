
#include <numbers>
#include <random>
#include <map>
#include <vector>
#include <unordered_map>

#include "IndexSin.h"
#include "VarDbl.h"
#include "Stat.h"

#ifndef __FFT_h__
#define __FFT_h__
namespace var_dbl 
{

class FFT {
public:
    enum SinSource {
        IndexedSin,
        LibrarySin
    };
    constexpr static const std::array<SinSource, 2> sSinSource{IndexedSin, LibrarySin};


    constexpr static unsigned char MIN_ORDER = 2;
    constexpr static unsigned char MAX_ORDER = 20;
        // (1 << MAX_ORDER) == 2 * PI
    VarDbl sin(long long freq, unsigned char order) const;
        // sin(Math.pi * freq / (1 << order))
    VarDbl cos(long long freq, unsigned char order) const;
        // cos(Math.pi * freq / (1 << order))

    static std::vector<size_t> bitReversedIndices(unsigned char order);

    static const IndexSin isin;

    std::vector<VarDbl> transform(const std::vector<VarDbl>& sData, bool forward) const;

    template<typename T> requires std::floating_point<T> || std::integral<T>
    std::vector<VarDbl> transform(const std::vector<T>& sData, bool forward) const;

    FFT(SinSource sinSource) { _sinType = sinSource; }
private:
    SinSource _sinType;
};


const IndexSin FFT::isin(FFT::MAX_ORDER);


inline std::vector<size_t> FFT::bitReversedIndices(unsigned char order) 
{
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

inline VarDbl FFT::sin(long long freq, unsigned char order) const
{
    std::ostringstream oss;
    switch (_sinType) {
    case IndexedSin:
        if (order < MIN_ORDER) {
            oss << "The order " << order << " < " << MIN_ORDER << " for fft.sin()";
            throw std::invalid_argument(oss.str());
        }
        if (MAX_ORDER < order) {
            oss << "The order " << order << " > " << MAX_ORDER << " for fft.sin()";
            throw std::invalid_argument(oss.str());
        }
        return isin.sin(freq *(1L <<(MAX_ORDER - order + 1)));
    case LibrarySin:            
        return std::sin(std::numbers::pi *freq /(1L << (order - 1)));
    default:
        oss << "Unknown SinSource " <<  _sinType << " for fft.sin()";
        throw std::invalid_argument(oss.str());
    }
}

inline VarDbl FFT::cos(long long freq, unsigned char order) const
{
    std::ostringstream oss;
    switch (_sinType) {
    case IndexedSin:
        if (order < MIN_ORDER) {
            oss << "The order " << order << " < " << MIN_ORDER << " for fft.cos()";
            throw std::invalid_argument(oss.str());
        }
        if (MAX_ORDER < order) {
            oss << "The order " << order << " > " << MAX_ORDER << " for fft.cos()";
            throw std::invalid_argument(oss.str());
        }
        return isin.cos(freq *(1L <<(MAX_ORDER - order + 1)));
    case LibrarySin:            
        return std::cos(std::numbers::pi *freq /(1L << (order - 1)));
    default:
        oss << "Unknown SinSource " <<  _sinType << " for fft.sin()";
        throw std::invalid_argument(oss.str());
    }
}

/*
    * 1-dimentional Fast Fourier Transformation (FFT)
    * 
    * @param sData     an array of size (2<<order), with each datum contains (real, image)
    * 
    * @return      an array of size (2<<order), with each datum contains (real, image)
    */
inline std::vector<VarDbl> FFT::transform(const std::vector<VarDbl>& sData, bool forward) const
{
    int order = MIN_ORDER;
    for (; order <= MAX_ORDER; ++order) {
        if ((2 << order) == sData.size()) {
            break;
        }
    }
    if (order > MAX_ORDER) {
        std::ostringstream oss;
        oss << "Invalid input array size " << sData.size() << " for fft.transform()";
        throw std::invalid_argument(oss.str());
    }
    const unsigned size = 1 << order;

    std::vector<VarDbl> sRes(2 << order);
    
    const std::vector<size_t> sIndex = bitReversedIndices(order);
    for (int i = 0; i < sIndex.size(); i++) {
        const int j = sIndex[i];
        sRes[(i << 1)] = sData[j << 1];
        sRes[(i << 1) + 1] = sData[(j << 1) + 1];
    }

    for (int i = 0; i < (sIndex.size() - 1); i += 2 ) {
        const VarDbl rt = sRes[(i << 1)], it = sRes[(i << 1) + 1];
        sRes[(i << 1)] += sRes[(i << 1) + 2];
        sRes[(i << 1) + 1] += sRes[(i << 1) + 3];
        sRes[(i << 1) + 2] = rt - sRes[(i << 1) + 2];
        sRes[(i << 1) + 3] = it - sRes[(i << 1) + 3];
    }

    for (unsigned o = 2, k = 4; o <= order; ++o, k <<= 1) {
        for (long j = 0; j < (k >> 1); j++) {
            const VarDbl vcos = cos( j, o );
            const VarDbl vsin = sin( forward? j : -j, o );
            for (int i = 0; i < sIndex.size(); i += k ) {
                const int idx0 = (i + j) << 1;
                const int idx1 = idx0 + k;
                const VarDbl r1 = sRes[idx1];
                const VarDbl i1 = sRes[idx1 + 1];
        
                const VarDbl rd = r1 * vcos - i1 * vsin;
                const VarDbl id = i1 * vcos + r1 * vsin;
        
                sRes[idx1] = sRes[idx0] - rd;
                sRes[idx1 + 1] = sRes[idx0 + 1] - id;
                sRes[idx0] += rd;
                sRes[idx0 + 1] += id;
            }   // for( i
        }
    }
    
    if (!forward) {
        for (int i = 0; i < (sIndex.size() << 1); i ++ ) {
            sRes[i] *= 1.0/(1L << order);
        }
    }
    return sRes;
}

template<typename T> requires std::floating_point<T> || std::integral<T>
inline std::vector<VarDbl> FFT::transform(const std::vector<T>& sData, bool forward) const
{
    return transform(std::vector<VarDbl>{sData}, forward);
}



} // namespace var_dbl
#endif  // __FFT_h__