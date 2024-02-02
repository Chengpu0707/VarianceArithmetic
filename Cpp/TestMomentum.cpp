#include "Momentum.h"
#include "Test.h"
#include "ulp.h"

using namespace var_dbl;

int main() 
{
    constexpr static const auto _momentum = Momentum<200, 5>();
    test::assertEquals(  1 * 0.99999947194737793, _momentum.factor(0));
    test::assertEquals(  1 * 0.99998569318184616, _momentum.factor(2));
    test::assertEquals(  3 * 0.99987013368935596, _momentum.factor(4));
    test::assertEquals( 15 * 0.99928863789470035, _momentum.factor(6));
    test::assertEquals(105 * 0.99719860134891214, _momentum.factor(8));
    test::assertEquals(945 * 0.99135593485973217, _momentum.factor(10)); 

    test::assertEquals(0, _momentum.factor(1));
    test::assertEquals(0, _momentum.factor(3));
    test::assertEquals(0, _momentum.factor(5));
    test::assertEquals(0, _momentum.factor(7));
    test::assertEquals(0, _momentum.factor(9));

    std::cout << "All Momentum tests are successful";
    return 0;
}