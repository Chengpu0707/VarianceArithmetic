#include "VarDbl.h"
#include "Test.h"
#include "ulp.h"

#include <functional>

using namespace var_dbl;


void assertEquals(VarDbl var, double value, double uncertainty,
                  double deltaValue=0, double deltaUncertainty=0) 
{
    if (deltaValue == 0)
        deltaValue = std::max(ulp(var.value()), ulp(value));
    test::assertAlmostEqual(var.value(), value, deltaValue);
    if (deltaUncertainty == 0)
        deltaUncertainty = std::max(ulp(var.uncertainty()), ulp(uncertainty));
    test::assertAlmostEqual(var.uncertainty(), uncertainty, deltaUncertainty);
}

void assertValueError(std::function<void()> func, std::string what ) 
{
    try {
        func();
        test::fail(what + std::string(" value overflow "));
    } catch(InitException ex) {
    } catch (std::exception ex) {
        test::fail(what + std::string(" catch not InitException but exception ") + ex.what());
    } catch (...) {
        test::fail(what + std::string(" catch unknown exception other than InitException"));
    }
}

void assertUncertaintyError(std::function<void()> func, std::string what ) 
{
    try {
        func();
        test::fail(what + std::string(" uncertainty overflow "));
    } catch(InitException ex) {
    } catch (std::exception ex) {
        test::fail(what + std::string(" catch not InitException but exception ") + ex.what());
    } catch (...) {
        test::fail(what + std::string(" catch unknown exception other than InitException"));
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
    long long ll = (1LL << 53); 
    assertEquals( VarDbl(ll), ll, 0 );
    ll = (1LL << 53) + 1;
    assertEquals( VarDbl(ll), ll, 0.5 );
    ll = (1LL << 54) + 1;
    assertEquals( VarDbl(ll), ll, 0.25 );
    ll = (1LL << 54) + 3;
    assertEquals( VarDbl(ll), ll, 0.25 );
}

void testInitFloat() {
    assertEquals( VarDbl(0.0f), 0, 0 );
    assertEquals( VarDbl(1.0f), 1, 0 );
    assertValueError([](){ VarDbl(std::numeric_limits<float>::infinity()); }, "VarDbl(float inf)");
}

void testInitDouble() {
    // use ulp 
    assertEquals( VarDbl(0.0), 0, 0 );
    assertEquals( VarDbl(-1.0), -1, VarDbl::ulp(std::numeric_limits<float>::epsilon()) );
    assertEquals( VarDbl(-2.0), -2, VarDbl::ulp(2 * std::numeric_limits<float>::epsilon()) );
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

    // InitException has precedence
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
    assertUncertaintyError([maxU](){ VarDbl(std::numeric_limits<double>::max(), maxU + ulp(maxU)); }, 
                         "VarDbl(double 1.79769e+308, double 1.34078e+154)");
    assertEquals( VarDbl(0, maxU), 0, maxU );
   
    const double minU = sqrt(ulp(std::numeric_limits<double>::min()));
    assertEquals( VarDbl(0, minU), 0, minU );
    assertEquals( VarDbl(0, minU/2), 0, 0 );
}

void testRepresentation() 
{
    const VarDbl v(1e-11, sqrt(2));
    test::assertTrue(v.to_string() == "1.000000e-11~1.414214e+00", v.to_string());
    std::ostringstream os;
    os.precision(17);
    os << v;
    test::assertEqual(os.str(), std::string("9.99999999999999939e-12~1.41421356237309515e+00"), os.str());
    std::istringstream is(os.str());
    VarDbl vr;
    is >> vr;
    test::assertEqual(vr.value(), v.value());
    test::assertEqual(vr.uncertainty(), v.uncertainty());
    os.str("");
    os << 1 << "+-" << 2;
    std::istringstream is2(os.str());
    try {
        is2 >> vr;
        test::fail("Illegal VarDbl stream accepted");
    } catch (std::invalid_argument ex) {
    } catch (...) {
        test::fail("Illegal VarDbl stream with invalid exception");
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
    assertEquals( v2 + 2.0, 3, 0 );
    assertEquals( 2.0 + v2, 3, 0 );
    v2 += 2.0;
    assertEquals( v2, 3, 0 );
}

void testSubFloat() 
{
    VarDbl v1(1, sqrt(2));
    assertEquals( v1 - 2.0, -1, sqrt(2) );
    assertEquals( 2.0 - v1, 1, sqrt(2) );
    v1 -= 2.0;
    assertEquals( v1, -1, sqrt(2) );

    VarDbl v2(1.0);
    assertEquals( v2 - 2.0, -1, 0 );
    assertEquals( 2.0 - v2, 1, 0 );
    v2 -= 2.0;
    assertEquals( v2, -1, 0 );
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

    assertEquals( VarDbl(-1.0) * VarDbl(2), -2, 0);
    assertEquals( VarDbl(-1.0) * 2, -2, 0);
    assertEquals( 2 * VarDbl(-1.0), -2, 0);

    assertEquals( VarDbl(-1.0) * VarDbl(2.0), -2,0);
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
    assertEquals(VarDbl(1) / VarDbl(1, 1e-3), 1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(VarDbl(-1) / VarDbl(1, 1e-3), -1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(VarDbl(1) / VarDbl(-1, 1e-3), -1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(VarDbl(-1) / VarDbl(-1, 1e-3), 1, 1e-3, 1e-6, 4.0e-9);
    assertEquals(VarDbl(2) / VarDbl(1, 1e-3), 2, 0.002, 2e-6, 6.4e-9);
    assertEquals(VarDbl(0.5) / VarDbl(1, 1e-3), 0.5, 0.0005, 5e-7, 4.0e-9);

    assertEquals(VarDbl(0, 1e-3) / VarDbl(1, 1e-3), 0, 1e-3, 0, 1.5e-9);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(1, 1e-3), 1, 1e-3*std::sqrt(2), 1e-6, 1.2e-9);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(1, 1e-3), -1, 1e-3*std::sqrt(2), 1e-6, 1.2e-9);
    assertEquals(VarDbl(1, 1e-3) / VarDbl(-1, 1e-3), -1, 1e-3*std::sqrt(2), 1e-6, 1.2e-9);
    assertEquals(VarDbl(-1, 1e-3) / VarDbl(-1, 1e-3), 1, 1e-3*std::sqrt(2), 1e-6, 1.2e-9);
    assertEquals(VarDbl(2, 1e-3) / VarDbl(1, 1e-3), 2, 1e-3*std::sqrt(5), 2e-6, 5.0e-9);
    assertEquals(VarDbl(0.5, 1e-3) / VarDbl(1, 1e-3), 0.5, 0.0005*std::sqrt(5), 5e-7, 1.2e-9);
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


void testCompareVarDbl()
{
    const VarDbl v1(1.000, 0.002);
    const VarDbl v2(1.001, 0.001);
    test::assertTrue(v1 == v2);
    test::assertFalse(v1 != v2);
    test::assertFalse(v1 < v2);
    test::assertTrue(v1 <= v2);
    test::assertFalse(v1 > v2);
    test::assertTrue(v1 >= v2);

    const VarDbl v3(1.002, 0.001);
    test::assertFalse(v1 == v3);
    test::assertTrue(v1 != v3);
    test::assertTrue(v1 < v3);
    test::assertTrue(v1 <= v3);
    test::assertFalse(v1 > v3);
    test::assertFalse(v1 >= v3);
    test::assertFalse(v3 < v1);
    test::assertFalse(v3 <= v1);
    test::assertTrue(v3 > v1);
    test::assertTrue(v3 >= v1);
}

void testCompareFloat()
{
    const VarDbl v1(1.000, 0.002);
    const float v2 = 1.001;
    test::assertTrue(v1 == v2);
    test::assertFalse(v1 != v2);
    test::assertFalse(v1 < v2);
    test::assertTrue(v1 <= v2);
    test::assertFalse(v1 > v2);
    test::assertTrue(v1 >= v2);

    test::assertTrue(v2 == v1);
    test::assertFalse(v2 != v1);
    test::assertFalse(v2 < v1);
    test::assertTrue(v2 <= v1);
    test::assertFalse(v2 > v1);
    test::assertTrue(v2 >= v1);

    const float v3 = 1.002;
    test::assertFalse(v1 == v3);
    test::assertTrue(v1 != v3);
    test::assertTrue(v1 < v3);
    test::assertTrue(v1 <= v3);
    test::assertFalse(v1 > v3);
    test::assertFalse(v1 >= v3);
    test::assertFalse(v3 < v1);
    test::assertFalse(v3 <= v1);
    test::assertTrue(v3 > v1);
    test::assertTrue(v3 >= v1);
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

    testVarDblByVarDblOne();
    testVarDblByFloatOne();
    testFloatByVarDblOne();
    testVarDblByVarDblZero();
    testVarDblByVarDblTwo();

    testCompareVarDbl();
    testCompareFloat();

    std::cout << "All VarDbl tests are successful";
    return 0;
}