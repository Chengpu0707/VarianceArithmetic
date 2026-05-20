#ifndef __Matrix_h__
#define __Matrix_h__

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#if __cplusplus >= 201703L
#include <filesystem>
#endif

#include "Stat.h"
#include "VarDbl.h"

namespace var_dbl
{

/**
 * Square matrix of VarDbl elements with determinant + adjugate computed via
 * bottom-up Laplace expansion for values (Formula 6.1) and direct evaluation
 * of Formula 6.9 for variances (independent-cell assumption).
 *
 * The static `runMatrixAnalysis` driver sweeps random + Hilbert matrices
 * across (size, noise) cells and writes MatrixCondition_*.txt and
 * AdjMatrix_*.txt under the given output directory. Cells are processed in
 * parallel via `std::async`; results are drained in submission order so the
 * output stays deterministic.
 */
class Matrix {
    const size_t        _size;
    std::vector<VarDbl> _sValue;
    VarDbl              _determ;
    std::vector<VarDbl> _sAdjugate;   // empty until calc()

    Matrix(size_t size, const std::vector<VarDbl>& sValue)
        : _size(size), _sValue(sValue) {}

public:
    explicit Matrix(size_t size)
        : _size(size), _sValue(size * size, VarDbl(0)) {}

    size_t        size() const                              { return _size; }
    const VarDbl& get(size_t row, size_t col) const         { return _sValue[row * _size + col]; }
    void          set(size_t row, size_t col, const VarDbl& v) {
        _sValue[row * _size + col] = v;
        _sAdjugate.clear();
    }

    VarDbl determ();
    Matrix adjugate();
    Matrix multiply(const Matrix& other) const;

    // Static factories
    static std::vector<std::vector<int> >
        createIntMatrix(size_t size, std::mt19937& rng, int range);
    static Matrix asMatrix    (const std::vector<std::vector<int> >& sInt);
    static Matrix addNoise    (const std::vector<std::vector<int> >& sInt, double dev, std::mt19937& rng);
    static Matrix hilbertMatrix(size_t size, double dev, std::mt19937& rng);
    static Matrix randomMatrix(size_t size, double dev);

    // Analysis driver
    static std::vector<double> defaultNoises();
    static void runMatrixAnalysis(size_t minSize, size_t maxSize, size_t targetCells,
                                  const std::vector<double>& noises,
                                  const std::string& outDir);

private:
    void calc();
};


// ----------------------------------------------------------------------------
//  Bit helpers
// ----------------------------------------------------------------------------

#if defined(__GNUC__) || defined(__clang__)
inline size_t _ctzll(unsigned long long x)      { return __builtin_ctzll(x); }
inline size_t _popcountll(unsigned long long x) { return __builtin_popcountll(x); }
#else
inline size_t _ctzll(unsigned long long x) {
    size_t c = 0; while (!(x & 1ULL)) { x >>= 1; ++c; } return c;
}
inline size_t _popcountll(unsigned long long x) {
    size_t c = 0; for (; x; x >>= 1) c += size_t(x & 1ULL); return c;
}
#endif

/** Recursive σ²-submatrix permanent — sum over m! bijections from set-bits of
 *  rowsMask to set-bits of colsMask of the product of m cell variances. */
inline double _permanent(unsigned long long rowsMask, unsigned long long colsMask,
                          const std::vector<double>& cellVar, size_t size) {
    if (rowsMask == 0) return 1.0;
    const size_t row = _ctzll(rowsMask);
    const unsigned long long newRows = rowsMask & ~(1ULL << row);
    const size_t rowOff = row * size;
    double sum = 0.0;
    unsigned long long rest = colsMask;
    while (rest) {
        const unsigned long long bit = rest & ~(rest - 1);
        const size_t col = _ctzll(bit);
        rest ^= bit;
        const double sig = cellVar[rowOff + col];
        if (sig != 0.0) sum += sig * _permanent(newRows, colsMask ^ bit, cellVar, size);
    }
    return sum;
}


// ----------------------------------------------------------------------------
//  calc(): determinant + adjugate via Formula (6.1) value and (6.9) variance
// ----------------------------------------------------------------------------

inline VarDbl Matrix::determ()   { if (_sAdjugate.empty()) calc(); return _determ; }
inline Matrix Matrix::adjugate() { if (_sAdjugate.empty()) calc(); return Matrix(_size, _sAdjugate); }

inline void Matrix::calc() {
    _sAdjugate.clear();
    if (_size == 0) { _determ = VarDbl(0); return; }

    using Mask = unsigned long long;
    const size_t nm = (_size >= 8 * sizeof(Mask)) ? size_t(-1) : (size_t(1) << _size);

    // (1) Group masks by popcount and build the mask→index lookup.
    std::vector<size_t> countAt(_size + 1, 0);
    for (Mask m = 0; m < nm; ++m) ++countAt[_popcountll(m)];
    std::vector<std::vector<Mask> > masksAtCount(_size + 1);
    for (size_t k = 0; k <= _size; ++k) masksAtCount[k].reserve(countAt[k]);
    for (Mask m = 0; m < nm; ++m) masksAtCount[_popcountll(m)].push_back(m);
    std::vector<size_t> maskIdx(nm);
    for (size_t k = 0; k <= _size; ++k)
        for (size_t i = 0; i < masksAtCount[k].size(); ++i)
            maskIdx[masksAtCount[k][i]] = i;

    // (2) Cell means and variances as primitive arrays.
    std::vector<double> cellVal(_size * _size), cellVar(_size * _size);
    for (size_t idx = 0; idx < _size * _size; ++idx) {
        cellVal[idx] = _sValue[idx].value();
        cellVar[idx] = _sValue[idx].variance();
    }

    // (3) Bottom-up Laplace on doubles: sub[k][rmIdx * cnt + cmIdx] = submatrix det.
    std::vector<std::vector<double> > sub(_size + 1);
    for (size_t k = 0; k <= _size; ++k) {
        const size_t c = masksAtCount[k].size();
        sub[k].assign(c * c, 0.0);
    }
    sub[0][0] = 1.0;
    for (size_t k = 1; k <= _size; ++k) {
        const std::vector<Mask>& masks = masksAtCount[k];
        const size_t cnt   = masks.size();
        const size_t cntm1 = masksAtCount[k - 1].size();
        const std::vector<double>& subKm1 = sub[k - 1];
        std::vector<double>& subK = sub[k];
        for (size_t ai = 0; ai < cnt; ++ai) {
            const Mask rm = masks[ai];
            const size_t firstRow  = _ctzll(rm);
            const size_t rmRestIdx = maskIdx[rm & ~(Mask(1) << firstRow)];
            const size_t rowOff    = firstRow * _size;
            for (size_t bi = 0; bi < cnt; ++bi) {
                const Mask cm = masks[bi];
                double det = 0.0;
                size_t posJ = 0;
                Mask rest = cm;
                while (rest) {
                    const Mask bit = rest & ~(rest - 1);
                    const size_t j = _ctzll(bit);
                    const double term = cellVal[rowOff + j]
                                       * subKm1[rmRestIdx * cntm1 + maskIdx[cm ^ bit]];
                    det += ((posJ++ & 1) == 0) ? term : -term;
                    rest ^= bit;
                }
                subK[ai * cnt + bi] = det;
            }
        }
    }

    // Formula (6.9) over a (containerRowMask × containerColMask) submatrix.
    auto formula69 = [&](Mask cRow, Mask cCol) {
        const size_t nC = _popcountll(cRow);
        double variance = 0.0;
        for (size_t m = 1; m <= nC; ++m) {
            const std::vector<Mask>& masks = masksAtCount[m];
            const size_t compLevel = nC - m;
            const size_t cnt = masksAtCount[compLevel].size();
            const std::vector<double>& subC = sub[compLevel];
            for (Mask posR : masks) {
                if ((posR & ~cRow) != 0) continue;
                const size_t rmIdx = maskIdx[cRow ^ posR];
                for (Mask posC : masks) {
                    if ((posC & ~cCol) != 0) continue;
                    const double sv = subC[rmIdx * cnt + maskIdx[cCol ^ posC]];
                    const double sv2 = sv * sv;
                    if (sv2 == 0.0) continue;
                    const double perm = _permanent(posR, posC, cellVar, _size);
                    if (perm > 0.0) variance += sv2 * perm;
                }
            }
        }
        return variance;
    };

    // (4) Determinant value (Formula 6.1) + variance via Formula (6.9).
    const Mask fullMask  = nm - 1;
    const size_t fullIdx = maskIdx[fullMask];
    const size_t cntN    = masksAtCount[_size].size();
    _determ = VarDbl(sub[_size][fullIdx * cntN + fullIdx],
                     std::sqrt(formula69(fullMask, fullMask)));

    // (5) Adjugate: value = ±(N−1)×(N−1) sub-det; variance from Formula (6.9).
    _sAdjugate.assign(_size * _size, VarDbl(0));
    const std::vector<double>& subKm1Final = sub[_size - 1];
    const size_t cntKm1 = masksAtCount[_size - 1].size();
    for (size_t i = 0; i < _size; ++i) {
        for (size_t j = 0; j < _size; ++j) {
            const Mask rmAdj = fullMask & ~(Mask(1) << j);
            const Mask cmAdj = fullMask & ~(Mask(1) << i);
            const double minor = subKm1Final[maskIdx[rmAdj] * cntKm1 + maskIdx[cmAdj]];
            const double adjVal = ((i + j) & 1) == 0 ? minor : -minor;
            _sAdjugate[i * _size + j] = VarDbl(adjVal, std::sqrt(formula69(rmAdj, cmAdj)));
        }
    }
}


// ----------------------------------------------------------------------------
//  Multiplication and static factories
// ----------------------------------------------------------------------------

inline Matrix Matrix::multiply(const Matrix& other) const {
    if (other._size != _size)
        throw std::invalid_argument("Matrix::multiply size mismatch");
    std::vector<VarDbl> out(_size * _size, VarDbl(0));
    for (size_t i = 0; i < _size; ++i) {
        const size_t rowOff = i * _size;
        for (size_t j = 0; j < _size; ++j) {
            double accVal = 0.0, accVar = 0.0;
            for (size_t k = 0; k < _size; ++k) {
                const VarDbl& a = _sValue[rowOff + k];
                const VarDbl& b = other._sValue[k * _size + j];
                const double av = a.value(),    bv = b.value();
                const double aV = a.variance(), bV = b.variance();
                accVal += av * bv;
                accVar += aV * bv * bv + bV * av * av + aV * bV;
            }
            out[i * _size + j] = VarDbl(accVal, std::sqrt(accVar));
        }
    }
    return Matrix(_size, out);
}

inline std::vector<std::vector<int> >
Matrix::createIntMatrix(size_t size, std::mt19937& rng, int range) {
    std::uniform_int_distribution<int> dist(-range, +range);
    std::vector<std::vector<int> > out(size, std::vector<int>(size));
    for (size_t r = 0; r < size; ++r)
        for (size_t c = 0; c < size; ++c) out[r][c] = dist(rng);
    return out;
}

inline Matrix Matrix::asMatrix(const std::vector<std::vector<int> >& sInt) {
    const size_t n = sInt.size();
    std::vector<VarDbl> sv(n * n, VarDbl(0));
    for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < n; ++c)
            sv[r * n + c] = VarDbl(static_cast<long long>(sInt[r][c]));
    return Matrix(n, sv);
}

inline Matrix Matrix::addNoise(const std::vector<std::vector<int> >& sInt,
                                double dev, std::mt19937& rng) {
    const size_t n = sInt.size();
    std::normal_distribution<double> gauss(0.0, dev);
    std::vector<VarDbl> sv(n * n, VarDbl(0));
    for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < n; ++c)
            sv[r * n + c] = VarDbl(sInt[r][c] + gauss(rng), dev);
    return Matrix(n, sv);
}

inline Matrix Matrix::hilbertMatrix(size_t size, double dev, std::mt19937& rng) {
    std::normal_distribution<double> gauss(0.0, dev > 0 ? dev : 1.0);
    std::vector<VarDbl> sv(size * size, VarDbl(0));
    for (size_t r = 0; r < size; ++r) {
        for (size_t c = 0; c < size; ++c) {
            double v = 1.0 / static_cast<double>(r + c + 1);
            const double unc = VarDbl::ulp(v);
            if (dev > 0) v += gauss(rng);
            sv[r * size + c] = VarDbl(v, std::sqrt(unc * unc + dev * dev));
        }
    }
    return Matrix(size, sv);
}

inline Matrix Matrix::randomMatrix(size_t size, double dev) {
    std::mt19937 rng(std::random_device{}());
    const double sigma = 256.0 * dev;
    std::uniform_int_distribution<int>  ui(-128, 128);
    std::normal_distribution<double>   gn(0.0, sigma);
    std::vector<VarDbl> sv(size * size, VarDbl(0));
    for (size_t idx = 0; idx < size * size; ++idx) {
        sv[idx] = VarDbl(ui(rng) + gn(rng), sigma);
    }
    return Matrix(size, sv);
}


// ----------------------------------------------------------------------------
//  runMatrixAnalysis driver (port of Python/runMatrixAnalysis.py)
// ----------------------------------------------------------------------------

namespace _matrix_analysis {

constexpr int    ELEMENT_RANGE = 1 << 8;
const     double DEV_SCALE     = ELEMENT_RANGE / std::sqrt(3.0);

// Histogram for normalized Adj error: 31 bins from -3 to +3 in 0.2 steps.
constexpr int    HISTO_DIVIDS = 5;
constexpr double HISTO_DEVS   = 3.0;
constexpr int    HISTO_HALF   = 15;       // round(HISTO_DIVIDS * HISTO_DEVS)
constexpr int    HISTO_BINS   = 2 * HISTO_HALF + 1;

inline void accumHisto(double v, int* histo) {
    if (v < -HISTO_DEVS || v > HISTO_DEVS) return;
    histo[HISTO_HALF + static_cast<int>(std::lround(v * HISTO_DIVIDS))]++;
}

inline std::string condHeader() {
    return "Size\tType\tNoise\tCondition Number"
           "\tDeterminant Value\tDeterminant Uncertainty"
           "\tDeterminant Precision\tRun Time\n";
}

inline std::string adjHeader() {
    static const char* CAT[] = {"Adj", "Fwd", "Rnd"};
    static const char* MET[] = {"Unc", "Val", "Norm"};
    static const char* FLD[] = {"Deviation", "Mean", "Minimum", "Maximum", "Count", "Loss"};
    std::ostringstream o;
    o << "Type\tNoise\tSize\tCount";
    for (int c = 0; c < 3; ++c)
        for (int m = 0; m < 3; ++m)
            for (int f = 0; f < 6; ++f)
                o << '\t' << CAT[c] << ' ' << MET[m] << ' ' << FLD[f];
    // Adj-norm histogram bin centers: -3.0, -2.8, …, 2.8, 3.0
    for (int i = -HISTO_HALF; i <= HISTO_HALF; ++i)
        o << '\t' << (static_cast<double>(i) / HISTO_DIVIDS);
    o << '\n';
    return o.str();
}

/** 1-norm condition: ||A||_1 * ||adj(A)||_1 / |det|. Order of magnitude matches
 *  numpy's 2-norm cond without pulling in an SVD dependency. */
inline double computeCond(const Matrix& m, const Matrix& adj, const VarDbl& det) {
    const double dv = det.value();
    if (dv == 0.0) return std::numeric_limits<double>::infinity();
    const size_t n = m.size();
    double normM = 0.0, normA = 0.0;
    for (size_t c = 0; c < n; ++c) {
        double sM = 0.0, sA = 0.0;
        for (size_t r = 0; r < n; ++r) {
            sM += std::abs(m.get(r, c).value());
            sA += std::abs(adj.get(r, c).value());
        }
        if (sM > normM) normM = sM;
        if (sA > normA) normA = sA;
    }
    return normM * normA / std::abs(dv);
}

inline void appendCondRow(std::ostringstream& f, size_t size, const char* type,
                           double noise, Matrix m) {
    using clk = std::chrono::high_resolution_clock;
    const auto t0 = clk::now();
    Matrix adj = m.adjugate();
    VarDbl det = m.determ();
    double cond;
    try { cond = computeCond(m, adj, det); }
    catch (...) { cond = std::numeric_limits<double>::quiet_NaN(); }
    const double val = det.value();
    const double unc = std::sqrt(det.variance());
    const double prec = (val != 0.0)
                        ? unc / std::abs(val)
                        : std::numeric_limits<double>::infinity();
    const double rt = std::chrono::duration<double>(clk::now() - t0).count();
    f << std::setprecision(std::numeric_limits<double>::max_digits10)
      << size << '\t' << type << '\t' << noise << '\t' << cond << '\t'
      << val << '\t' << unc << '\t' << prec << '\t' << rt << '\n';
}

inline void feedStats(Stat<double> stats[3], int losses[3],
                      double unc, double val, double normUnc) {
    stats[0].add(unc);
    stats[1].add(val);
    if (normUnc > 0.0) stats[2].add(val / normUnc);
    else               losses[2]++;
}

/** The Adj normalization uses an UNCONDITIONAL predictor — Formula (6.9)
 *  applied to a (sInt means, dev uncertainty) matrix — so the recorded
 *  `Adj Norm` reflects "how well does Formula (6.9) predict the actual
 *  noise-induced spread of adj(M)?". Fwd and Rnd still use the
 *  sample-conditional uncertainty propagation through M·adj. */
inline void accumulate(const std::vector<std::vector<int> >& sInt, double dev,
                       Matrix mNoisy, const Matrix& adjPrecise,
                       Stat<double> stats[3][3], int losses[3][3],
                       int* adjHisto) {
    const size_t n = mNoisy.size();

    // Unconditional Adj predictor.
    Matrix preciseWithDev(n);
    for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < n; ++c)
            preciseWithDev.set(r, c, VarDbl(static_cast<double>(sInt[r][c]), dev));
    Matrix adjPredicted = preciseWithDev.adjugate();

    Matrix adjNoisy = mNoisy.adjugate();
    VarDbl detNoisy = mNoisy.determ();
    const double detVal = detNoisy.value();
    const double detVar = detNoisy.variance();
    Matrix ssId  = mNoisy.multiply(adjNoisy);
    Matrix ssIdL = adjNoisy.multiply(mNoisy);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const VarDbl& adjCell = adjNoisy.get(i, j);
            const VarDbl& preCell = adjPrecise.get(i, j);
            const VarDbl& idCell  = ssId.get(i, j);
            const VarDbl& idLCell = ssIdL.get(i, j);

            // Adj: unconditional predictor → ratio ≈ 1 if Formula (6.9) is faithful.
            const double adjPredUnc = std::sqrt(adjPredicted.get(i, j).variance());
            const double adjErr     = adjCell.value() - preCell.value();
            feedStats(stats[0], losses[0], adjPredUnc, adjErr, adjPredUnc);
            if (adjPredUnc > 0.0) accumHisto(adjErr / adjPredUnc, adjHisto);

            // Fwd: (M·adj − det·I)[i,j] — rounding error of an identity; the
            // sample-conditional propagated unc is the right denominator.
            const double idVal = idCell.value();
            const double idUnc = std::sqrt(idCell.variance());
            const double fwdUnc = (i == j) ? std::sqrt(idUnc * idUnc + detVar) : idUnc;
            feedStats(stats[1], losses[1], fwdUnc,
                      (i == j) ? idVal - detVal : idVal, fwdUnc);

            // Rnd: (M·adj − adj·M)[i,j] — same reasoning as Fwd.
            const double idLVal = idLCell.value();
            const double idLUnc = std::sqrt(idLCell.variance());
            const double rndUnc = std::sqrt(idUnc * idUnc + idLUnc * idLUnc);
            feedStats(stats[2], losses[2], rndUnc, idVal - idLVal, rndUnc);
        }
    }
}

inline std::pair<std::string, std::string>
runCell(size_t size, double noise, size_t count) {
    const double dev = DEV_SCALE * noise;
    static thread_local std::mt19937 rng(static_cast<unsigned>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
        ^ std::hash<std::thread::id>{}(std::this_thread::get_id())));

    std::ostringstream cond;
    try { appendCondRow(cond, size, "Hilbert", noise, Matrix::hilbertMatrix(size, dev, rng)); }
    catch (...) {}

    Stat<double> stats[3][3];
    int          losses[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    int          adjHisto[HISTO_BINS] = {0};
    size_t actual = 0;
    for (size_t r = 0; r < count; ++r) {
        try {
            const auto sInt = Matrix::createIntMatrix(size, rng, ELEMENT_RANGE);
            Matrix adjPrecise = Matrix::asMatrix(sInt).adjugate();
            Matrix mNoisy = Matrix::addNoise(sInt, dev, rng);
            appendCondRow(cond, size, "Random", noise, mNoisy);
            accumulate(sInt, dev, std::move(mNoisy), adjPrecise, stats, losses, adjHisto);
            ++actual;
        } catch (...) {}
    }

    std::ostringstream adj;
    adj << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "Gaussian\t" << noise << '\t' << size << '\t' << actual;
    for (int c = 0; c < 3; ++c) {
        for (int m = 0; m < 3; ++m) {
            Stat<double>& s = stats[c][m];
            adj << '\t' << s.std()  << '\t' << s.mean()
                << '\t' << s.min()  << '\t' << s.max()
                << '\t' << s.count() << '\t' << losses[c][m];
        }
    }
    // 31 bin probabilities for the normalized adj error.
    long long total = 0;
    for (int i = 0; i < HISTO_BINS; ++i) total += adjHisto[i];
    for (int i = 0; i < HISTO_BINS; ++i) {
        adj << '\t';
        if (total > 0) adj << (static_cast<double>(adjHisto[i]) / static_cast<double>(total));
    }
    adj << '\n';
    return {cond.str(), adj.str()};
}

inline std::set<std::string> readDoneCells(const std::string& adjPath) {
    std::set<std::string> done;
    std::ifstream prev(adjPath);
    if (!prev) return done;
    std::string line;
    std::getline(prev, line);    // header
    while (std::getline(prev, line)) {
        std::istringstream iss(line);
        std::string tok, noiseStr, sizeStr;
        if (!std::getline(iss, tok, '\t'))      continue;
        if (!std::getline(iss, noiseStr, '\t')) continue;
        if (!std::getline(iss, sizeStr, '\t'))  continue;
        if (!noiseStr.empty() && !sizeStr.empty())
            done.insert(sizeStr + "|" + noiseStr);
    }
    return done;
}

} // namespace _matrix_analysis

inline std::vector<double> Matrix::defaultNoises() {
    std::vector<double> v(20);
    v[0] = 0.0;
    for (int i = 1; i < 20; ++i) v[i] = std::pow(10.0, i - 18);
    return v;
}

inline void Matrix::runMatrixAnalysis(size_t minSize, size_t maxSize, size_t targetCells,
                                       const std::vector<double>& noises,
                                       const std::string& outDir) {
    if (minSize >= maxSize)
        throw std::invalid_argument("runMatrixAnalysis: empty size range");
#if __cplusplus >= 201703L
    std::filesystem::create_directories(outDir);
#endif
    const std::string condPath = outDir + "/MatrixCondition_"
                              + std::to_string(minSize) + "_" + std::to_string(maxSize) + ".txt";
    const std::string adjPath  = outDir + "/AdjMatrix_"
                              + std::to_string(minSize) + "_" + std::to_string(maxSize) + ".txt";

    const std::set<std::string> done = _matrix_analysis::readDoneCells(adjPath);
    if (!done.empty())
        std::cout << "Skipping " << done.size() << " existing cells from " << adjPath << '\n';

    const bool condExisted = std::ifstream(condPath).good();
    const bool adjExisted  = std::ifstream(adjPath ).good();
    std::ofstream fc(condPath, condExisted ? std::ios::app : std::ios::trunc);
    std::ofstream fa(adjPath,  adjExisted  ? std::ios::app : std::ios::trunc);
    if (!condExisted) fc << _matrix_analysis::condHeader();
    if (!adjExisted)  fa << _matrix_analysis::adjHeader();

    std::vector<std::future<std::pair<std::string, std::string> > > futures;
    for (size_t size = minSize; size < maxSize; ++size) {
        const size_t count = std::max<size_t>(1, targetCells / (size * size));
        for (double noise : noises) {
            std::ostringstream key;
            key << std::setprecision(std::numeric_limits<double>::max_digits10);
            key << size << "|" << noise;
            if (done.count(key.str())) continue;
            futures.emplace_back(std::async(std::launch::async,
                [size, noise, count]() { return _matrix_analysis::runCell(size, noise, count); }));
        }
    }
    for (auto& f : futures) {
        try {
            auto r = f.get();
            fc << r.first;  fc.flush();
            fa << r.second; fa.flush();
        } catch (const std::exception& ex) {
            std::cerr << "cell failed: " << ex.what() << '\n';
        }
    }
}

}  // namespace var_dbl

#endif
