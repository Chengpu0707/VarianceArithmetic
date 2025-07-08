#include "TestTaylor.h"

#include <functional>

using namespace var_dbl;

void testVarDblByVarDblOne() 
{
    assertEquals(VarDbl(0) / VarDbl(1), 0, 0);
    assertEquals(VarDbl(1) / VarDbl(1),   1, 0);
    assertEquals(VarDbl(-1) / VarDbl(1), -1, 0);
    assertEquals(VarDbl(1) / VarDbl(-1), -1, 0);
    assertEquals(VarDbl(-1) / VarDbl(-1), 1, 0);
    assertEquals(VarDbl(2) / VarDbl(1),   2, 0);

    assertEquals(VarDbl(1.0) / VarDbl(1), 1, 0);
    assertEquals(VarDbl(2.0) / VarDbl(1), 2, 0);
    assertEquals(VarDbl(0.5) / VarDbl(1), 0.5, 0);

    assertEquals(VarDbl(1.0) / VarDbl(1.0), 1, 0, 0, 0);
    assertEquals(VarDbl(2.0) / VarDbl(1.0), 2, 0, 0, 0);
    assertEquals(VarDbl(0.5) / VarDbl(1.0), 0.5, 0, 0, 0);

    assertEquals(VarDbl(0, 1e-3) / VarDbl(1), 0, 1e-3);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(1), 1, 1e-3);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(1), -1, 1e-3);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(-1), -1, 1e-3);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(-1), 1, 1e-3);
    assertEquals(VarDbl(2, 1e-3) / VarDbl(1), 2, 1e-3);
    assertEquals(VarDbl(0.5, 1e-3) / VarDbl(1), 0.5, 1e-3);

    assertEquals(VarDbl(0) / VarDbl(1, 1e-3), 0, 0);
    assertEquals(VarDbl(1) / VarDbl(1, 1e-3), 1, 1e-3, 1e-6, 4e-9);
    assertEquals(VarDbl(-1) / VarDbl(1, 1e-3), -1, 1e-3, 1e-6, 4e-9);
    assertEquals(VarDbl(1) / VarDbl(-1, 1e-3), -1, 1e-3, 1e-6, 4e-9);
    assertEquals(VarDbl(-1) / VarDbl(-1, 1e-3), 1, 1e-3, 1e-6, 4e-9);
    assertEquals(VarDbl(2) / VarDbl(1, 1e-3), 2, 0.002, 2e-6, 8e-9);
    assertEquals(VarDbl(0.5) / VarDbl(1, 1e-3), 0.5, 0.0005, 5e-7, 4e-9);

    assertEquals(VarDbl(0, 1e-3) / VarDbl(1, 1e-3), 0, 1e-3, 0, 2e-9);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(1, 1e-3), 1, 1e-3*std::sqrt(2), 1e-6, 2e-9);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(1, 1e-3), -1, 1e-3*std::sqrt(2), 1e-6, 2e-9);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(-1, 1e-3), -1, 1e-3*std::sqrt(2), 1e-6, 2e-9);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(-1, 1e-3), 1, 1e-3*std::sqrt(2), 1e-6, 2e-9);
    assertEquals(VarDbl(2, 1e-3) / VarDbl(1, 1e-3), 2, 1e-3*std::sqrt(5), 2e-6, 6e-9);
    assertEquals(VarDbl(0.5, 1e-3) / VarDbl(1, 1e-3), 0.5, 0.0005*std::sqrt(5), 5e-7, 2e-9);
}

void testVarDblByFloatOne() 
{
    assertEquals(VarDbl(0)/ 1.0, 0, 0);
    assertEquals(VarDbl(1)/ 1.0,   1, 0);
    assertEquals(VarDbl(-1)/ 1.0, -1, 0);
    assertEquals(VarDbl(1)/ -1.0, -1, 0);
    assertEquals(VarDbl(-1)/ -1.0, 1, 0);
    assertEquals(VarDbl(2)/ 1.0,   2, 0);

    assertEquals(VarDbl(0.5)/ 1.0, 0.5, 0);
    assertEquals(VarDbl(1.0)/ 1.0, 1.0, 0);
    assertEquals(VarDbl(2.0)/ 1.0, 2.0, 0);

    assertEquals(VarDbl(0, 1e-3)/ 1.0, 0, 1e-3);
    assertEquals(VarDbl(1, 1e-3)/ 1.0, 1, 1e-3);
    assertEquals(VarDbl(-1, 1e-3)/ 1.0, -1, 1e-3);
    assertEquals(VarDbl(1, 1e-3)/ -1.0, -1, 1e-3);
    assertEquals(VarDbl(-1, 1e-3)/ -1.0, 1, 1e-3);
    assertEquals(VarDbl(2, 1e-3)/ 1.0, 2, 1e-3);
    assertEquals(VarDbl(0.5, 1e-3)/ 1.0, 0.5, 1e-3);
}

void testFloatByVarDblOne() 
{
    assertEquals(0 / VarDbl(1), 0, 0);
    assertEquals(1 / VarDbl(1),   1, 0);
    assertEquals(-1 / VarDbl(1), -1, 0);
    assertEquals(1 / VarDbl(-1), -1, 0);
    assertEquals(-1 / VarDbl(-1), 1, 0);
    assertEquals(2 / VarDbl(1),   2, VarDbl::ulp(2.));

    assertEquals(0.5 / VarDbl(1), 0.5, 0);

    assertEquals(0 / VarDbl(1, 1e-3), 0, 0);
    assertEquals(1 / VarDbl(1, 1e-3), 1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(-1 / VarDbl(1, 1e-3), -1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(1 / VarDbl(-1, 1e-3), -1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(-1 / VarDbl(-1, 1e-3), 1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(2 / VarDbl(1, 1e-3), 2, 0.002, 2e-6, 6.4e-9);
    assertEquals(0.5 / VarDbl(1, 1e-3), 0.5, 0.0005, 5e-7, 4.0e-9);
}

void testVarDblByVarDblZero() 
{
    try {
        VarDbl(0) / VarDbl(0);
        test::fail("VarDbl(0) / VarDbl(0)");
    } catch (InitException ex) {
    } catch (std::exception ex) {
        test::fail(ex.what());
    }

    try {
        VarDbl(0) / VarDbl(0, 1);
        test::fail("VarDbl(0) / VarDbl(0, 1)");
    } catch (InitException) {
    } catch (std::exception ex) {
        test::fail(ex.what());
    }

    try {
        VarDbl(0) / VarDbl(0.1, 1);
        test::fail("VarDbl(0) / VarDbl(0.1, 1)");
    } catch (NotMonotonicException ex) {
    } catch (NotReliableException ex) {
    } catch (std::exception ex) {
        test::fail(ex.what());
    }
}

void testVarDblByVarDblTwo() 
{
    assertEquals(VarDbl(0) / VarDbl(2), 0, 0);
    assertEquals(VarDbl(1) / VarDbl(2),   0.5, 0);
    assertEquals(VarDbl(-1) / VarDbl(2), -0.5, 0);
    assertEquals(VarDbl(1) / VarDbl(-2), -0.5, 0);
    assertEquals(VarDbl(-1) / VarDbl(-2), 0.5, 0);
    assertEquals(VarDbl(2) / VarDbl(2),     1, 0);
    assertEquals(VarDbl(0.5) / VarDbl(2), 0.25, 0);

    assertEquals(VarDbl(0, 1e-3) / VarDbl(2), 0, 5e-4);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(2), 0.5, 5e-4);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(2), -0.5, 5e-4);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(-2), -0.5, 5e-4);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(-2), 0.5, 5e-4);
    assertEquals(VarDbl(2, 1e-3) / VarDbl(2), 1, 5e-4);
    assertEquals(VarDbl(0.5, 1e-3) / VarDbl(2), 0.25, 5e-4);

    assertEquals(VarDbl(0) / VarDbl(2, 1e-3), 0, 0);
    assertEquals(VarDbl(1) / VarDbl(2, 1e-3), 0.5, 2.5e-4, 1e-6, 4.0e-9);
    assertEquals(VarDbl(-1) / VarDbl(2, 1e-3), -0.5, 2.5e-4, 1e-6, 4.0e-9);
    assertEquals(VarDbl(1) / VarDbl(-2, 1e-3), -0.5, 2.5e-4, 1e-6, 4.0e-9);
    assertEquals(VarDbl(-1) / VarDbl(-2, 1e-3), 0.5, 2.5e-4, 1e-6, 4.0e-9);
    assertEquals(VarDbl(2) / VarDbl(2, 1e-3), 1, 5e-4, 2e-6, 6.4e-9);
    assertEquals(VarDbl(0.5) / VarDbl(2, 1e-3), 0.25, 1.25e-4, 5e-7, 4.0e-9);

    assertEquals(VarDbl(0, 1e-3) / VarDbl(2, 1e-3), 0, 5e-4, 0, 1.5e-9);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(2, 1e-3), 0.5, 2.5e-4*std::sqrt(5), 1e-6, 1.2e-9);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(2, 1e-3), -0.5, 2.5e-4*std::sqrt(5), 1e-6, 1.2e-9);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(-2, 1e-3), -0.5, 2.5e-4*std::sqrt(5), 1e-6, 1.2e-9);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(-2, 1e-3), 0.5, 2.5e-4*std::sqrt(5), 1e-6, 1.2e-9);
    assertEquals(VarDbl(2, 1e-3) / VarDbl(2, 1e-3), 1, 1e-3*std::sqrt(0.5), 2e-6, 5.0e-9);
    assertEquals(VarDbl(0.5, 1e-3) / VarDbl(2, 1e-3), 0.25, 2.5e-4*std::sqrt(4.25), 5e-7, 1.2e-9);
}


int main() 
{
    testVarDblByVarDblOne();
    testVarDblByFloatOne();
    testFloatByVarDblOne();
    testVarDblByVarDblZero();
    testVarDblByVarDblTwo();
}