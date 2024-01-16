#include "ValDbl.h"

#include <cassert>
#include <iostream>


void testInitInt() {
    const VarDbl v0(0), v1(-1);
    assert(v0.value() == 0);
    assert(v0.uncertainty() == 0);
    assert(v1.value() == -1);
    assert(v1.uncertainty() == 0);

    // when long long looses resolution
    const long long lU = -(1LL << 53); 
    const VarDbl vU(lU), vO(lU - 1);
    assert(vU.value() == lU);
    assert(vU.uncertainty() == 0);
    assert(vO.value() == (double) (lU - 1));
    assert(vO.uncertainty() == 2 / sqrt(3));
}

void testInitDouble() {
    // use ulp 
    const VarDbl v0(0.0), v1(-1.0), v1p(-1, 0);
    assert(v0.value() == 0);
    assert(v0.uncertainty() == 0);
    assert(v1.value() == -1);
    assert(v1.uncertainty() == std::numeric_limits<double>::epsilon() /sqrt(3));
    assert(v1p.value() == -1);
    assert(v1p.uncertainty() == 0);
}

void testValueError(double value, std::string what ) 
{
    try {
        const VarDbl v(value);
        assert(false);
    } catch(ValueError ex) {
        assert(ex.what == what);
    } catch (std::exception ex) {
        assert(false);
    } catch (...) {
        assert(false);
    }
}

void testValueError(double value, double uncertainty, std::string what ) 
{
    try {
        const VarDbl v(value, uncertainty);
        assert(false);
    } catch(ValueError ex) {
        assert(ex.what == what);    
    } catch (std::exception ex) {
        assert(false);
    } catch (...) {
        assert(false);
    }
}

void testUncertaintyError(double value, double uncertainty, std::string what ) 
{
    try {
        const VarDbl v(value, uncertainty);
        assert(false);
    } catch(UncertaintyError ex) {
        assert(ex.what == what);    
    } catch (std::exception ex) {
        assert(false);
    } catch (...) {
        assert(false);
    }
}

void testInitNaN() {
    testValueError(1.0 / 0.0, "VarDbl(double value)");
    testValueError(std::numeric_limits<double>::infinity(), 
                    "VarDbl(double value)");
    testValueError(std::numeric_limits<double>::quiet_NaN(), "VarDbl(double value)");

    // ValueError has precedence
    testValueError(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                    "VarDbl(double value, double uncertainty)");
    testValueError(std::numeric_limits<double>::infinity(), 0,
                    "VarDbl(double value, double uncertainty)");

    testUncertaintyError(0, std::numeric_limits<double>::quiet_NaN(), 
                        "VarDbl(double value, double uncertainty)");
    // variance calculation overflow
    testUncertaintyError(0, std::numeric_limits<double>::max(), 
                        "VarDbl(double value, double uncertainty)");
}

void testInit() {
    testInitInt();
    testInitDouble();
    testInitNaN();
}

int main() {
    testInit();

    std::cout << "All VarDbl tests are successful";
    return 0;
}