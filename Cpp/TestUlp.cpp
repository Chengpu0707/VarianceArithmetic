/*
Unit tests for ulp.h: verifies units-in-the-last-place computation for integer
and floating-point inputs and confirms boundary behaviour at the precision
limits of double.

Cross-language probe (invoked with argument "Test"): prints the VarDbl
uncertainty deviation for representative integer boundaries, to verify ULP
consistency between C++ and the Java UlpProbeJava reference.
*/
#include <cassert>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>

#include "ulp.h"
#include "VarDbl.h"

using namespace var_dbl;


void testUlp()
{
    assert(0 == var_dbl::ulp(0));
    assert(0 == var_dbl::ulp(1));
    assert(0 == var_dbl::ulp(-1));
    long long max = 1LL << std::numeric_limits<double>::digits;
    assert(0 == var_dbl::ulp(max - 1));
    assert(0 == var_dbl::ulp(max));
    assert(1 == var_dbl::ulp(max + 1));
    long long u = (max << 1) + max + 3;
    assert(-0.5 == var_dbl::ulp(u));
    assert(-1 == u - ((long long) ((double) u)));

    assert(std::numeric_limits<double>::denorm_min() == var_dbl::ulp(0.));
    assert(std::numeric_limits<double>::denorm_min() == var_dbl::ulp(std::numeric_limits<double>::min()));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(1.));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(-1.));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(std::sqrt(2)));
    assert(std::numeric_limits<double>::epsilon()*2 == var_dbl::ulp(2.));
    assert(std::numeric_limits<double>::epsilon()*0.5 == var_dbl::ulp(0.5));

    assert(0 == var_dbl::ulp(0., 2));
    assert(0 == var_dbl::ulp(1., 2));
    assert(0 == var_dbl::ulp(-1., 2));
    assert(0 == var_dbl::ulp(-1., 2));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(std::sqrt(2), 2));
}


// Cross-language probe: prints VarDbl uncertainty deviation for representative
// integer values at and beyond the 2^53 boundary, where double-precision loses
// integer exactness. Output is compared line-by-line with Java's UlpProbeJava
// (Java/bin/UlpProbeJava.class) to verify ULP consistency across languages.
void dumpUlp()
{
    long long vals[] = {(1LL<<53)+1LL, (1LL<<53)+3LL, (1LL<<53)+5LL,
                        (1LL<<53)+7LL, (1LL<<53)+15LL,
                        (1LL<<54)+1LL, (1LL<<55)+5LL,
                        7LL*(1LL<<50), 127LL*(1LL<<50)};
    for (long long v : vals) {
        VarDbl x(v);
        printf("C++  ulp(0x%016llX) = %g\n", v, std::sqrt(x.variance()));
    }
}


int main(int argc, char* argv[])
{
    if (argc == 1) {
        testUlp();
        std::cout << "All ulp tests are successful";
    } else if ((argc == 2) && (std::string(argv[1]) == "Test"))
        dumpUlp();
    return 0;
}
