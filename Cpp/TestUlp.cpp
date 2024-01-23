#include <cassert>

#include "ulp.h"

int main() 
{
    assert(4.9406564584124654e-324 == var_dbl::ulp(0));
    assert(4.9406564584124654e-324 == var_dbl::ulp(std::numeric_limits<double>::min()));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(1));
    assert(std::numeric_limits<double>::epsilon() == var_dbl::ulp(-1));
}