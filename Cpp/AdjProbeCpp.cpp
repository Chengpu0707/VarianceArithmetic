#include <cstdio>
#include <cmath>
#include "Matrix.h"
using namespace var_dbl;
int main() {
    std::vector<std::vector<int> > sInt(8, std::vector<int>(8, 0));
    long long seed = 1;
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) {
            seed = seed * 6364136223846793005LL + 1442695040888963407LL;
            sInt[i][j] = (int)((seed >> 32) % 257) - 128;
        }
    Matrix m = Matrix::asMatrix(sInt);
    Matrix adj = m.adjugate();
    printf("C++ adj[0,0] value=%.10e unc=%.10e\n", adj.get(0,0).value(), std::sqrt(adj.get(0,0).variance()));
    Matrix prod = m.multiply(adj);
    printf("C++ (M·adj)[0,0] value=%.10e unc=%.10e\n", prod.get(0,0).value(), std::sqrt(prod.get(0,0).variance()));
    printf("C++ (M·adj)[0,1] value=%.10e unc=%.10e\n", prod.get(0,1).value(), std::sqrt(prod.get(0,1).variance()));
    Matrix m2 = Matrix::asMatrix(sInt);
    VarDbl det = m2.determ();
    printf("C++ det value=%.10e unc=%.10e\n", det.value(), std::sqrt(det.variance()));
    return 0;
}
