/*
Smoke test that simply includes every public header in the library to verify
they compile cleanly together and that there are no missing transitive includes
across all supported C++ standards.
*/
#include <cstdlib>
#include <iostream>

// test all headers
#include "ulp.h"
#include "FFT.h"
#include "IndexSin.h"
#include "Stat.h"
#include "Taylor.h"
#include "Test.h"
#include "TestFFT.h"
#include "VarDbl.h"

using namespace var_dbl;



int main()
{
    std::cout << "All header tests are successful";
    return 0;
}