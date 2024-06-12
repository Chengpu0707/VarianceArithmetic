#include <cassert>
#include <cmath>
#include <iostream>

#include "ulp.h"

int main() 
{
    assert(0 == var_dbl::ulp(0));
    assert(0 == var_dbl::ulp(1));
    assert(0 == var_dbl::ulp(-1));
    long long max = 1LL << std::numeric_limits<double>::digits;
    assert(0 == var_dbl::ulp(max - 1));
    assert(0 == var_dbl::ulp(max));
    assert(0.5 == var_dbl::ulp(max + 1));
    long long u = (max << 1) + max + 3;
    assert(-0.25 == var_dbl::ulp(u));
    assert(-1 == u - ((long long) ((double) u)));

    assert(std::numeric_limits<double>::denorm_min() == var_dbl::ulp(0.));
    assert(std::numeric_limits<double>::denorm_min() == var_dbl::ulp(std::numeric_limits<double>::min()));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(1.));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(-1.));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(std::sqrt(2)));

    assert(0 == var_dbl::ulp(0., 2));
    assert(0 == var_dbl::ulp(1., 2));
    assert(0 == var_dbl::ulp(-1., 2));
    assert(0 == var_dbl::ulp(-1., 2));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(std::sqrt(2), 2));

    std::cout << "All ulp tests are successful";
}