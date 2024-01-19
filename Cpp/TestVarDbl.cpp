#include "ValDbl.h"

using namespace var_dbl;

void testInitCopy() {
    const VarDbl vs(-1, 4);
    const VarDbl vc(vs);
    Test::assertEqual(vs.value(), vc.value());
    Test::assertEqual(vs.uncertainty(), vc.uncertainty());
}

void testInitInt() {
    const VarDbl v0, v1(-1);
    Test::assertEqual(v0.value(), 0.0);
    Test::assertEqual(v0.uncertainty(), 0.0);
    Test::assertEqual(v1.value(), -1.0);
    Test::assertEqual(v1.uncertainty(),  0.0);

    // when long long looses resolution
    const long long lU = -(1LL << 53); 
    const VarDbl vU(lU), vO(lU - 1);
    Test::assertEqual(vU.value(), (double) lU);
    Test::assertEqual(vU.uncertainty(), 0.0);
    Test::assertEqual(vO.value(), (double) (lU - 1));
    Test::assertEqual(vO.uncertainty(), 2 / sqrt(3));
}

void testInitFloat() {
    const VarDbl v0(0.0f), v1(-1.0f);
    Test::assertEqual(v0.value(), 0.0);
    Test::assertEqual(v0.uncertainty(), 8.0904004559294244e-46);
        // std::numeric_limits<float>::min() == 1.1754943508222875e-38
    Test::assertEqual(v1.value(), -1.0);
    Test::assertEqual(v1.uncertainty(), 
        std::numeric_limits<float>::epsilon() * VarDbl::DEVIATION_OF_LSB);

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
    const VarDbl v0(0.0), v1(-1.0), v1p(-1, 0);
    Test::assertEqual(v0.value(), 0.0);
    Test::assertEqual(v0.uncertainty(), 0.0);
    Test::assertEqual(v1.value(), -1.0);
    Test::assertEqual(v1.uncertainty(), 
        std::numeric_limits<double>::epsilon() * VarDbl::DEVIATION_OF_LSB);
    Test::assertEqual(v1p.value(), -1.0);
    Test::assertEqual(v1p.uncertainty(), 0.0);
}

void testInitUncertainty() {
    const VarDbl v0(1, 0), v1(-1, 1), v1p(1, -1);
    Test::assertEqual(v0.value(), 1.0);
    Test::assertEqual(v0.uncertainty(), 0.0);
    Test::assertEqual(v1.value(), -1.0);
    Test::assertEqual(v1.uncertainty(), 1.0);
    Test::assertEqual(v1p.value(), 1.0);
    Test::assertEqual(v1p.uncertainty(), 1.0);
}

void testValueError(double value, std::string what ) 
{
    try {
        const VarDbl v(value);
        Test::fail();
    } catch(ValueError ex) {
        Test::assertEqual(ex.what, what);
    } catch (std::exception ex) {
        Test::fail();
    } catch (...) {
        Test::fail();
    }
}

void testValueError(double value, double uncertainty, std::string what ) 
{
    try {
        const VarDbl v(value, uncertainty);
        Test::fail();
    } catch(ValueError ex) {
        Test::assertEqual(ex.what, what);    
    } catch (std::exception ex) {
        Test::fail();
    } catch (...) {
        Test::fail();
    }
}

void testUncertaintyError(double value, double uncertainty, std::string what ) 
{
    try {
        const VarDbl v(value, uncertainty);
        Test::fail();
    } catch(VarianceError ex) {
        Test::assertEqual(ex.what, what);    
    } catch (std::exception ex) {
        Test::fail();
    } catch (...) {
        Test::fail();
    }
}

void testInitException() {
    testValueError(1.0 / 0.0, "VarDbl(double inf)");
    testValueError(std::numeric_limits<double>::infinity(), 
                    "VarDbl(double inf)");
    testValueError(std::numeric_limits<double>::quiet_NaN(), "VarDbl(double nan)");

    // ValueError has precedence
    testValueError(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                    "VarDbl(double inf, double inf)");
    testValueError(std::numeric_limits<double>::infinity(), 0,
                    "VarDbl(double inf, double 0)");

    testUncertaintyError(0, std::numeric_limits<double>::quiet_NaN(), 
                        "VarDbl(double 0, double nan)");
    // variance calculation overflow
    testUncertaintyError(0, std::numeric_limits<double>::max(), 
                        "VarDbl(double 0, double 1.79769e+308)");
}

void testUncertaintyRange() {
    const double maxU = sqrt(std::numeric_limits<double>::max());
    const double minU = sqrt(Test::ulp(std::numeric_limits<double>::min()));

    testUncertaintyError(std::numeric_limits<double>::max(), maxU + Test::ulp(maxU), 
                         "VarDbl(double 1.79769e+308, double 1.34078e+154)");
    
    const VarDbl max(std::numeric_limits<double>::max(), maxU - Test::ulp(maxU));
    Test::assertEquals(max.value(), std::numeric_limits<double>::max());
    Test::assertEquals(max.uncertainty(), maxU);
    
    const VarDbl min(std::numeric_limits<double>::min(), minU);
    Test::assertEquals(min.value(), std::numeric_limits<double>::min());
    Test::assertEquals(min.uncertainty(), minU);
    
    const VarDbl zero(std::numeric_limits<double>::min(), minU * 0.5);
    Test::assertEquals(zero.value(), std::numeric_limits<double>::min());
    Test::assertEquals(zero.uncertainty(), 0);
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
    VarDbl v = v1 + v2;
    Test::assertEquals(0, v.value());
    Test::assertEquals(2, v.uncertainty());

    v1 += v2;
    Test::assertEquals(0, v1.value());
    Test::assertEquals(2, v1.uncertainty());
}

void testSubVarDbl() 
{
    VarDbl v1(sqrt(2), sqrt(2));
    VarDbl v2(-sqrt(2), sqrt(2));

    VarDbl v = v2 - v1;
    Test::assertEquals(-2*sqrt(2), v.value());
    Test::assertEquals(2, v.uncertainty());

    v2 -= v1;
    Test::assertEquals(-2*sqrt(2), v2.value());
    Test::assertEquals(2, v2.uncertainty());
}

void testAddInt() 
{
    VarDbl v1(1, sqrt(2));
    VarDbl v = v1 + 2;
    Test::assertEquals(3, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    v = 2 + v1;
    Test::assertEquals(3, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    v1 += 2;
    Test::assertEquals(3, v1.value());
    Test::assertEquals(sqrt(2), v1.uncertainty());
}

void testSubInt() 
{
    VarDbl v1(1, sqrt(2));
    VarDbl v = v1 - 2;
    Test::assertEquals(-1, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    v = 2 - v1;
    Test::assertEquals(1, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    v1 -= 2;
    Test::assertEquals(-1, v1.value());
    Test::assertEquals(sqrt(2), v1.uncertainty());
}

void testAddFloat() 
{
    VarDbl v1(1, sqrt(2));
    VarDbl v = v1 + 2.0;
    Test::assertEquals(3, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    v = 2.0 + v1;
    Test::assertEquals(3, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    VarDbl v2(1.0);
    v = v2 + 2.0;
    Test::assertEquals(3, v.value());
    Test::assertTrue(Test::ulp(3.5) < v.uncertainty() < Test::ulp(4.0));

    v = 2.0 + v2;
    Test::assertEquals(3, v.value());
    Test::assertTrue(Test::ulp(3.5) < v.uncertainty() < Test::ulp(4.0));

    v1 += 2.0;
    Test::assertEquals(3, v1.value());
    Test::assertEquals(sqrt(2), v1.uncertainty());

    v2 += 2.0;
    Test::assertEquals(3, v2.value());
    Test::assertTrue(Test::ulp(3.5) < v2.uncertainty() < Test::ulp(4.0));
}

void testSubFloat() 
{
    VarDbl v1(1, sqrt(2));
    VarDbl v = v1 - 2.0;
    Test::assertEquals(-1, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    v = 2.0 - v1;
    Test::assertEquals(1, v.value());
    Test::assertEquals(sqrt(2), v.uncertainty());

    VarDbl v2(1.0);
    v = v2 - 2.0;
    Test::assertEquals(-1, v.value());
    Test::assertTrue(Test::ulp(3.5) < v.uncertainty() < Test::ulp(4.0));

    v = 2.0 - v2;
    Test::assertEquals(1, v.value());
    Test::assertTrue(Test::ulp(3.5) < v.uncertainty() < Test::ulp(4.0));

    v1 -= 2.0;
    Test::assertEquals(-1, v1.value());
    Test::assertEquals(sqrt(2), v1.uncertainty());

    v2 -= 2.0;
    Test::assertEquals(-1, v2.value());
    Test::assertTrue(Test::ulp(3.5) < v2.uncertainty() < Test::ulp(4.0));
}

void testAddSubException()
{
    const double maxV = std::numeric_limits<double>::max();
    const double maxU = sqrt(maxV);
    try {
        VarDbl(maxV, maxU) + VarDbl(maxV, maxU);
        Test::fail("testAddSubException() value overflow");
    } catch (ValueError ex) { 
    } catch (...) {
        Test::fail("testAddSubException() catch wrong ValueError");
    }

    try {
        VarDbl(maxV, maxU) - VarDbl(maxV, maxU);
        Test::fail("testAddSubException() variance overflow");
    } catch (VarianceError ex) { 
    } catch (...) {
        Test::fail("testAddSubException() catch wrong VarianceError");
    }
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

    std::cout << "All VarDbl init tests are successful";
    return 0;
}