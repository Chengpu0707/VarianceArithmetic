#include <cassert>
#include <iostream>

#include "ulp.h"

int main() 
{
    assert(4.9406564584124654e-324 == var_dbl::ulp(0));
    assert(4.9406564584124654e-324 == var_dbl::ulp(std::numeric_limits<double>::min()));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(1));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(-1));
    for (int i = 0; i < 10; ++i)
        std::cout << "1/(1+" << i/10.0 << ")=" << var_dbl::ulp(1/(1+i/10.0)) << '\n';
    for (int i = 0; i < 10; ++i)
        std::cout << "1/(1-" << i/10.0 << ")=" << var_dbl::ulp(1/(1-i/10.0)) << '\n';
}