#include "FFT.h"
#include "Test.h"

using namespace var_dbl;


void testBitReversion()
{
    test::assertEquals( FFT::bitReversedIndices(2), std::vector{0, 2, 1, 3}, "FFT::bitReversedIndices(2)" );
    test::assertEquals( FFT::bitReversedIndices(3), std::vector{0, 4, 2, 6, 1, 5, 3, 7}, "FFT::bitReversedIndices(3)");
    test::assertEquals( FFT::bitReversedIndices(4), 
                        std::vector{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15}, "FFT::bitReversedIndices(4)" );

}

int main()
{
    testBitReversion();
}

