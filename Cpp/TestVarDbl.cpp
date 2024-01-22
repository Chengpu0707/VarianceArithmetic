#include "Test.h"
#include "ValDbl.h"

#include <functional>

using namespace var_dbl;


void assertEquals(VarDbl var, double value, double uncertainty) 
{
    Test::assertEquals(var.value(), value);
    Test::assertEquals(var.uncertainty(), uncertainty);
}

void assertValueError(std::function<void()> func, std::string what ) 
{
    try {
        func();
        Test::fail(what + std::string(" value overflow "));
    } catch(ValueError ex) {
    } catch (std::exception ex) {
        Test::fail(what + std::string(" catch not ValueError but exception ") + ex.what());
    } catch (...) {
        Test::fail(what + std::string(" catch unknown exception other than ValueError"));
    }
}

void assertUncertaintyError(std::function<void()> func, std::string what ) 
{
    try {
        func();
        Test::fail(what + std::string(" uncertainty overflow "));
    } catch(UncertaintyError ex) {
    } catch (std::exception ex) {
        Test::fail(what + std::string(" catch not UncertaintyError but exception ") + ex.what());
    } catch (...) {
        Test::fail(what + std::string(" catch unknown exception other than UncertaintyError"));
    }
}





void testInitCopy() {
    const VarDbl vs(-1, 4);
    const VarDbl vc(vs);
    assertEquals( vc, vs.value(), vs.uncertainty() );
}

void testInitInt() {
    assertEquals( VarDbl(), 0, 0 );
    assertEquals( VarDbl(-1), -1, 0 );

    // when long long looses resolution
    const long long ll = (1LL << 53); 
    assertEquals( VarDbl(ll), ll, 0 );
    assertEquals( VarDbl(ll + 1), ll + 1, 2 / sqrt(3) );
    assertEquals( VarDbl(ll << 1), ll << 1, 0 );
}

void testInitFloat() {
    assertEquals( VarDbl(0.0f), 0, 8.0904004559294244e-46 );
        // std::numeric_limits<float>::min() == 1.1754943508222875e-38
    assertEquals( VarDbl(1.0f), 1, std::numeric_limits<float>::epsilon() * VarDbl::DEVIATION_OF_LSB );

    try {
        const VarDbl v(std::numeric_limits<float>::infinity());
        Test::fail();
    } catch(ValueError ex) {
    } catch (std::exception ex) {
        Test::fail();
    } catch (...) {
        Test::fail();
    }

}

void testInitDouble() {
    // use ulp 
    assertEquals( VarDbl(0.0), 0, 0 );
    assertEquals( VarDbl(-1.0), -1, std::numeric_limits<double>::epsilon() * VarDbl::DEVIATION_OF_LSB );
    assertEquals( VarDbl(-2.0), -2, 2 * std::numeric_limits<double>::epsilon() * VarDbl::DEVIATION_OF_LSB );
}

void testInitUncertainty() {
    assertEquals( VarDbl(1, 0), 1, 0 );
    assertEquals( VarDbl(-1, 1), -1, 1 );
    assertEquals( VarDbl(1, -1), 1, 1 );
}

void testInitException() {
    assertValueError([](){ VarDbl(1.0 / 0.0); }, "VarDbl(1.0 / 0.0)");
    assertValueError([](){ VarDbl(std::numeric_limits<double>::infinity()); }, "VarDbl(double inf)");
    assertValueError([](){ VarDbl(std::numeric_limits<double>::quiet_NaN()); }, "VarDbl(double nan)");

    // ValueError has precedence
    assertValueError([](){ VarDbl(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()); },
                    "VarDbl(double inf, double inf)");
    assertValueError([](){ VarDbl(std::numeric_limits<double>::infinity(), 0); },
                    "VarDbl(double inf, double 0)");

    assertUncertaintyError([](){ VarDbl(0, std::numeric_limits<double>::quiet_NaN()); }, 
                        "VarDbl(double 0, double nan)");
    // variance calculation overflow
    assertUncertaintyError([](){ VarDbl(0, std::numeric_limits<double>::max()); }, 
                        "VarDbl(double 0, double overflow)");
}

void testUncertaintyRange() {
    const double maxU = sqrt(std::numeric_limits<double>::max());
    assertUncertaintyError([maxU](){ VarDbl(std::numeric_limits<double>::max(), maxU + Test::ulp(maxU)); }, 
                         "VarDbl(double 1.79769e+308, double 1.34078e+154)");
    assertEquals( VarDbl(0, maxU), 0, maxU );
   
    const double minU = sqrt(VarDbl::ulp(std::numeric_limits<double>::min()));
    assertEquals( VarDbl(0, minU), 0, minU );
    assertEquals( VarDbl(0, minU/2), 0, 0 );
}

void testRepresentation() 
{
    const VarDbl v(1e-11, sqrt(2));
    Test::assertTrue(v.to_string() == "1.000000e-11~1.414214e+00", v.to_string());
    std::ostringstream os;
    os.precision(17);
    os << v;
    Test::assertEqual(os.str(), std::string("9.99999999999999939e-12~1.41421356237309515e+00"), os.str());
    std::istringstream is(os.str());
    VarDbl vr;
    is >> vr;
    Test::assertEquals(vr.value(), v.value());
    Test::assertEquals(vr.uncertainty(), v.uncertainty());
    os.str("");
    os << 1 << "+-" << 2;
    std::istringstream is2(os.str());
    try {
        is2 >> vr;
        Test::fail("Illegal VarDbl stream accepted");
    } catch (std::invalid_argument ex) {
    } catch (...) {
        Test::fail("Illegal VarDbl stream with invalid exception");
    }
}

void testAddVarDbl() 
{
    VarDbl v1(sqrt(2), sqrt(2));
    VarDbl v2(-sqrt(2), sqrt(2));
    assertEquals( v1 + v2, 0, 2 );
    v1 += v2;
    assertEquals( v1, 0, 2 );
    
}

void testSubVarDbl() 
{
    VarDbl v1(sqrt(2), sqrt(2));
    VarDbl v2(-sqrt(2), sqrt(2));
    assertEquals( v1 - v2, 2*sqrt(2), 2 );
    v2 -= v1;
    assertEquals( v2, -2*sqrt(2), 2 );
}

void testAddInt() 
{
    VarDbl v1(1, sqrt(2));
    assertEquals( v1 + 2, 3, sqrt(2) );
    assertEquals( 2 + v1, 3, sqrt(2) );
    v1 += 2;
    assertEquals( v1, 3, sqrt(2) );
}

void testSubInt() 
{
    VarDbl v1(1, sqrt(2));
    assertEquals( v1 - 2, -1, sqrt(2) );
    assertEquals( 2 - v1, 1, sqrt(2) );
    v1 -= 2;
    assertEquals( v1, -1, sqrt(2) );
}

void testAddFloat() 
{
    VarDbl v1(1, sqrt(2));
    assertEquals( v1 + 2.0, 3, sqrt(2) );
    assertEquals( 2.0 + v1, 3, sqrt(2) );
    v1 += 2.0;
    assertEquals( v1, 3, sqrt(2) );

    VarDbl v2(1.0);
    VarDbl v = v2 + 2.0;
    assertEquals( v2 + 2.0, 3, VarDbl::ulp(1)*sqrt(5) );
    assertEquals( 2.0 + v2, 3, VarDbl::ulp(1)*sqrt(5) );
    v2 += 2.0;
    assertEquals( v2, 3, VarDbl::ulp(1)*sqrt(5) );
}

void testSubFloat() 
{
    VarDbl v1(1, sqrt(2));
    assertEquals( v1 - 2.0, -1, sqrt(2) );
    assertEquals( 2.0 - v1, 1, sqrt(2) );
    v1 -= 2.0;
    assertEquals( v1, -1, sqrt(2) );

    VarDbl v2(1.0);
    assertEquals( v2 - 2.0, -1, VarDbl::ulp(1)*sqrt(5) );
    assertEquals( 2.0 - v2, 1, VarDbl::ulp(1)*sqrt(5) );
    v2 -= 2.0;
    assertEquals( v2, -1, VarDbl::ulp(1)*sqrt(5) );
}

void testAddSubException()
{
    const double maxV = std::numeric_limits<double>::max();
    const double maxU = sqrt(maxV);
    assertValueError([maxV, maxU](){VarDbl(maxV, maxU) + VarDbl(maxV, maxU);}, "testAddSubException()");
    assertUncertaintyError([maxV, maxU](){VarDbl(maxV, maxU) - VarDbl(maxV, maxU);}, "testAddSubException()");
}

void testMultiplyZero() 
{
    assertEquals( VarDbl(0) * VarDbl(1), 0, 0);
    assertEquals( VarDbl(0) * 1, 0, 0);
    assertEquals(  1 * VarDbl(0), 0, 0);

    assertEquals( VarDbl(0) * VarDbl(1.0), 0, 0);
    assertEquals( VarDbl(0) * 1.0, 0, 0);
    assertEquals( 1.0 * VarDbl(0), 0, 0);

    assertEquals( VarDbl(0, 1e-3) * VarDbl(2.0), 0, 2e-3);
    assertEquals( VarDbl(0, 1e-3) * 2.0, 0, 2e-3);
    assertEquals( 2.0 * VarDbl(0, 1e-3), 0, 2e-3);

    assertEquals( VarDbl(0) * VarDbl(2.0, 1e-3), 0, 0);
    assertEquals( VarDbl(2, 1e-3) * 0.0, 0, 0);
    assertEquals( 0.0 * VarDbl(2, 1e-3), 0, 0);

    assertEquals( VarDbl(0, 1e-3) * VarDbl(2.0), 0, 2e-3);
    assertEquals( VarDbl(0, 1e-3) * 2.0, 0, 2e-3);
    assertEquals( 2.0 * VarDbl(0, 1e-3), 0, 2e-3);

    assertEquals( VarDbl(0, 1e-3) * VarDbl(2.0, 1e-2), 0, sqrt(4 + 1e-4)*1e-3);
}

void testMultiplyOne() 
{
    assertEquals( VarDbl(-1) * VarDbl(2), -2, 0);
    assertEquals( VarDbl(-1) * 2, -2, 0);
    assertEquals( 2 * VarDbl(-1), -2, 0);

    assertEquals( VarDbl(-1.0) * VarDbl(2), -2, VarDbl::ulp(2.0));
    assertEquals( VarDbl(-1.0) * 2, -2, VarDbl::ulp(2.0));
    assertEquals( 2 * VarDbl(-1.0), -2, VarDbl::ulp(2.0));

    assertEquals( VarDbl(-1.0) * VarDbl(2.0), -2, sqrt(2)*VarDbl::ulp(2.0));
    assertEquals( VarDbl(-1.0, 1e-3) * VarDbl(2.0), -2, 2e-3);
    assertEquals( VarDbl(-1.0) * VarDbl(2.0, 1e-3), -2, 1e-3);
    assertEquals( VarDbl(-1.0, 1e-3) * VarDbl(2.0, 1e-3), -2, sqrt(5 + 1e-6) * 1e-3);
}

void testMultiplyException() 
{
    const double maxV = sqrt(std::numeric_limits<double>::max()) * 1.001;
    const double maxU = sqrt(std::numeric_limits<double>::max()) * 0.5;
    assertValueError([maxV, maxU](){VarDbl(maxV, maxU) * VarDbl(maxV, maxU);}, "testMultiplyException()");
    assertUncertaintyError([maxU](){VarDbl(1, maxU) * VarDbl(1, maxU);}, "testMultiplyException()");
}





int main() 
{
    testInitCopy();
    testInitInt();
    testInitFloat();
    testInitDouble();
    testInitUncertainty();
    testInitException();
    testUncertaintyRange();

    testRepresentation();

    testAddVarDbl();
    testSubVarDbl();
    testAddInt();
    testSubInt();
    testAddFloat();
    testSubFloat();

    testMultiplyZero();
    testMultiplyOne();
    testMultiplyException();

    std::cout << "All VarDbl init tests are successful";
    return 0;
}