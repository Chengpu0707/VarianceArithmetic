#include "ValDbl.h"



void testInitInt() {
    const VarDbl v0(0), v1(-1);
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

void testInitDouble() {
    // use ulp 
    const VarDbl v0(0.0), v1(-1.0), v1p(-1, 0);
    Test::assertEqual(v0.value(), 0.0);
    Test::assertEqual(v0.uncertainty(), 0.0);
    Test::assertEqual(v1.value(), -1.0);
    Test::assertEqual(v1.uncertainty(), std::numeric_limits<double>::epsilon() /sqrt(3));
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
    } catch(UncertaintyError ex) {
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
    const VarDbl v(-sqrt(2), sqrt(2));
    Test::assertTrue(v.to_string() == "-1.414214~1.414214", v.to_string());
    std::ostringstream os;
    os.precision(17);
    os << v;
    Test::assertEqual(os.str(), std::string("-1.4142135623730951~1.4142135623730951"), os.str());
    std::istringstream is(os.str());
    VarDbl vr;
    is >> vr;
    Test::assertEquals(vr.value(), v.value(), Test::ulp(v.value()));
    Test::assertEquals(vr.uncertainty(), v.uncertainty(), Test::ulp(v.uncertainty()));
}

int main() {
    testInitInt();
    testInitDouble();
    testInitUncertainty();
    testInitException();
    testUncertaintyRange();
    testRepresentation();

    std::cout << "All VarDbl init tests are successful";
    return 0;
}