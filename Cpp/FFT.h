

#include <unordered_map>
#include <vector>



#ifndef __fft_h__
#define __fft_h__
namespace var_dbl 
{

class FFT {
public:
    constexpr static unsigned MAX_ORDER = 20;
        // (1 << MAX_ORDER) == 2 * PI

    virtual double sin(int freq, unsigned char order) const = 0;
        // sin(Math.pi * freq / (1 << order))

    virtual double cos(int freq, unsigned char order) const = 0;
        // cos(Math.pi * freq / (1 << order))

    static std::vector<size_t> bitReversedIndices(unsigned char order);

};


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


} // namespace var_dbl
#endif  // __fft_h__