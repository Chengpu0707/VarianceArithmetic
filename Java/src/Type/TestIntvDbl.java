package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;

import java.util.HashMap;
import java.util.Map;

import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;


public class TestIntvDbl {
    IntvDbl op;
    IntvDbl op2;
    IntvDbl res;
    final double TOLERANCE = 2E-16;
    final double TOLERANCE_DIV = 5E-16;

    @Test 
    public void testSimple() {
        try {
            final double in = Math.sqrt(Double.MIN_NORMAL);
            multiply(in, 0, in, 0, false);
        } catch (ValueException | UncertaintyException | TypeException e) {
            fail();
        }
    }

    @Test
    public void testHash() {
        final Map<IntvDbl, Integer> sMap = new HashMap<IntvDbl, Integer>();
       try {
            final IntvDbl one = new IntvDbl(1, 0);
            final IntvDbl two = new IntvDbl(1, 2);
            final IntvDbl three = new IntvDbl(0, 1);
            sMap.put(one, 1);
            sMap.put(two, 2);
            sMap.put(three, 3);
            assertEquals((Integer) 1, sMap.get(one));
            assertEquals((Integer) 2, sMap.get(two));
            assertEquals((Integer) 3, sMap.get(three));
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    private void testToString(String str, double value, double range) {
        try {
            op = new IntvDbl(value, range);
            assertEquals(str, op.toString());
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testToString() {
        testToString("0", 0, 0);
        testToString("-1@100", -1, 100);
        testToString("-1@1.0e+03", -1, 1000);    
        testToString("1.000e-03@1.0e-03", 0.001, 0.001);
        testToString("-8.988e+307@9.0e+307", -Double.MAX_VALUE / 2, Double.MAX_VALUE / 2);
        testToString("2.225e-308@4.9e-324", Double.MIN_NORMAL, Double.MIN_VALUE);
        testToString("2.225e-308@1.8e+308", Double.MIN_NORMAL, Double.MAX_VALUE);        
    }


    private void testClone() {
        res = op.clone();
        assertTrue(op.equals(res));
        assertTrue(res.equals(op));
        assertTrue(op != res);
        assertFalse(op == res);

        IntvDbl clone = res.clone();
        assertTrue(op.equals(clone));
        assertTrue(clone.equals(op));
    }
    
    private void testNegate() {
        res = op.clone();
        res.negate();
        assertEquals(-op.value(), res.value(), TOLERANCE);
        assertEquals(op.uncertainty(), res.uncertainty(), TOLERANCE);
    }
    
    private void testShift( int bits ) {
        try {
            res = op.clone();
            res.shift(bits);
            if (bits >= 0) {
                assertEquals(op.value() * (1L << bits), res.value(), TOLERANCE);
                assertEquals(op.uncertainty() * (1L << bits), res.uncertainty(), TOLERANCE);
            } else {
                assertEquals(op.value() / (1L << -bits), res.value(), TOLERANCE);
                assertEquals(op.uncertainty() / (1L << -bits), res.uncertainty(), TOLERANCE);
            }
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testShift() {
        testShift(0);
        testShift(10);
        testShift(-10);
    }

    private void testOffset( double offset ) {
        try {
            res = op.clone();
            res.add(offset);
            assertEquals(op.value() + offset, res.value(), TOLERANCE);
            assertEquals(op.uncertainty(), res.uncertainty(), TOLERANCE);
        } catch (ValueException e) {
            fail(e.getMessage());
        }
    }
    private void testOffset() {
        testOffset(0);
        testOffset(10);
        testOffset(-10);
        testOffset(0.1);
        testOffset(-0.1);
        testOffset(1000);
        testOffset(-1000);
        testOffset(0.001);
        testOffset(-0.001);
    }

    private void testScale( double fold ) {
        try {
            res = op.clone();
            res.multiply(fold);
            assertEquals(op.value() * fold, res.value(), TOLERANCE);
            assertEquals(op.uncertainty() * Math.abs(fold), res.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testScale() {
        testScale(0);
        testScale(10);
        testScale(-10);
        testScale(0.1);
        testScale(-0.1);
        testScale(1000);
        testScale(-1000);
        testScale(0.001);
        testScale(-0.001);
    }

    private void testPowerAllRange( double exponent ) {
        try {
            res = op.clone();
            res.power(exponent);
            final double div = op.uncertainty() / 8;
            final double tolerance = (0 <= exponent)? TOLERANCE : TOLERANCE_DIV;
            if (div == 0) {
                assertEquals(Math.pow(op.value(), exponent), res.value(), tolerance);
            } else {
                if (((op.value() - op.uncertainty()) <= 0) && 
                    (0 <= (op.value() - op.uncertainty()))) {
                    assertTrue((res.value() - res.uncertainty()) <= 0);
                    assertTrue(0 <= (res.value() + res.uncertainty()));
                }
                for (int i = -8; i <= 8; ++i) {
                    final double d = Math.pow(op.value() + i*div, exponent);
                    if ((d + tolerance) < (res.value() - res.uncertainty()))
                        fail(String.format("(%e in %s)^%e: value %e < %s by %e", op.value() + i*div, op.toString(), exponent, 
                             d, res.toString(), (res.value() - res.uncertainty()) - d));
                    if ((res.value() + res.uncertainty()) < (d - tolerance)) {
                        fail(String.format("(%e in %s)^%e: value %e > %s by %e", op.value() + i*div, op.toString(), exponent, 
                             d, res.toString(), d - (res.value() + res.uncertainty())));
                    }
                }
            }
        } catch (ValueException e) {
            fail(e.getMessage());
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testPowerValueException( double exponent ) {
        try {
            res = op.clone();
            res.power(exponent);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testPowerUncertaintyException( double exponent ) {
        try {
            res = op.clone();
            res.power(exponent);
            fail();
        } catch (ValueException e) {
            fail(e.getMessage());
        } catch (UncertaintyException e) {
        }
    }
    private void testPower( double exponent, boolean raiseValueException ) {
        if ((Math.floor(exponent) != Math.ceil(exponent)) && 
            ((op.value() - op.uncertainty()) < 0)) {
            testPowerUncertaintyException( exponent );
            return;
        }
        if ((exponent < 0) && 
            ((op.value() - op.uncertainty()) <= 0) && 
            (0 <= (op.value() + op.uncertainty()))) {
            testPowerUncertaintyException( exponent );
            return;
        }
        if (exponent == 0) {
            try {
                res = op.clone();
                res.power(exponent);
                assertEquals(1, res.value(), 0);
                assertEquals(0, res.uncertainty(), 0);
            } catch (ValueException | UncertaintyException e) {
                fail(e.getMessage());
            }
            return;
        } 
        if (raiseValueException) {
            testPowerValueException(exponent);
        } else {
            testPowerAllRange( exponent );
        }
    }
    private void testPower(final double[] sFiniteExponent, final double[] sValueExceptionExponent) {
        if (sFiniteExponent != null) {
            for (double exp: sFiniteExponent) {
                testPower( exp, false );
            }
        }

        if (sValueExceptionExponent != null) {
            for (double exp: sValueExceptionExponent) {
                testPower( exp, true );
            }
        }
    }
    private void testPower(double value, double range, double exponent, boolean raiseValueException) {
        try {
            op = new IntvDbl(value, range);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
        testPower(exponent, raiseValueException);
    }

    private void testSingle(double value, double range, 
                            final double[] sFiniteExponent, final double[] sValueExceptionExponent) {
        try {
            op = new IntvDbl(value, range);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
        assertEquals(value, op.value(), TOLERANCE);
        assertEquals(range, op.uncertainty(), TOLERANCE);
        testClone();
        testNegate();
        testShift();
        testOffset();
        testScale();
        testPower(sFiniteExponent, sValueExceptionExponent);
    }
    private void testSingle(double value, double range) {
        testSingle(value, range, new double[]{0, 1, -1, 2, -2, 0.5, -0.5}, null);
    }

    @Test
    public void test0() {
        testSingle(0, 0);
        testSingle(0, 1);
        testSingle(0, 10);
        testSingle(0, 0.1);
        testSingle(0, 1000);
        testSingle(0, 0.001);
        testSingle( 0, Double.MIN_VALUE);
        testSingle( 0, Double.MIN_NORMAL);
        testSingle( 0, Math.sqrt(Double.MAX_VALUE));

        testPower(0, Double.MAX_VALUE, 0.5, false);
        testPower(0, Double.MAX_VALUE, -0.5, false);
        testPower(0, Double.MAX_VALUE, 2, true);
        testPower(0, Double.MAX_VALUE, -2, true);
    }

    @Test
    public void test1() {
        testSingle(1, 0);
        testSingle(-1, 0);
        testSingle(1, 1);
        testSingle(-1, 1);
        testSingle(1, 10);
        testSingle(-1, 10);
        testSingle(1, 0.1);
        testSingle(-1, 0.1);
        testSingle(1, 1000);
        testSingle(-1, 1000);
        testSingle(1, 0.001);
        testSingle(-1, 0.001);
    }

    @Test
    public void test4() {
        testSingle(4, 0);
        testSingle(-4, 0);
        testSingle(4, 1);
        testSingle(-4, 1);
        testSingle(4, 10);
        testSingle(-4, 10);
        testSingle(4, 0.1);
        testSingle(-4, 0.1);
        testSingle(4, 1000);
        testSingle(-4, 1000);
        testSingle(4, 0.001);
        testSingle(-4, 0.001);
    }

    @Test
    public void test4th() {
        testSingle(0.25, 0);
        testSingle(-0.25, 0);
        testSingle(0.25, 1);
        testSingle(-0.25, 1);
        testSingle(0.25, 10);
        testSingle(-0.25, 10);
        testSingle(0.25, 0.1);
        testSingle(-0.25, 0.1);
        testSingle(0.25, 1000);
        testSingle(-0.25, 1000);
        testSingle(0.25, 0.001);
        testSingle(-0.25, 0.001);
    }

    @Test
    public void testMinNormal() {
        final double[] sFiniteExponent = new double[] {0,1,-1,2,0.5,-0,5};
        final double[] sValueExceptionExponent = new double[]{-2};
        testSingle( Double.MIN_NORMAL, 0, sFiniteExponent, sValueExceptionExponent);
        testSingle(-Double.MIN_NORMAL, 0, sFiniteExponent, sValueExceptionExponent);
        testSingle( Double.MIN_NORMAL, Double.MIN_VALUE, sFiniteExponent, sValueExceptionExponent);
        testSingle(-Double.MIN_NORMAL, Double.MIN_VALUE, sFiniteExponent, sValueExceptionExponent);
        testSingle( Double.MIN_NORMAL, Double.MIN_NORMAL, sFiniteExponent, sValueExceptionExponent);
        testSingle(-Double.MIN_NORMAL, Double.MIN_NORMAL, sFiniteExponent, sValueExceptionExponent);
    }

    @Test
    public void testMinValue() {
        final double[] sFiniteExponent = new double[] {0,1,2,0.5,-0,5};
        final double[] sValueExceptionExponent = new double[]{-1, -2};
        testSingle( Double.MIN_VALUE, 0, sFiniteExponent, sValueExceptionExponent);
        testSingle(-Double.MIN_VALUE, 0, sFiniteExponent, sValueExceptionExponent);
        testSingle( Double.MIN_VALUE, Double.MIN_VALUE, sFiniteExponent, sValueExceptionExponent);
        testSingle(-Double.MIN_VALUE, Double.MIN_VALUE, sFiniteExponent, sValueExceptionExponent);
        testSingle( Double.MIN_VALUE, Double.MIN_NORMAL, sFiniteExponent, sValueExceptionExponent);
        testSingle(-Double.MIN_VALUE, Double.MIN_NORMAL, sFiniteExponent, sValueExceptionExponent);
    }

    @Test
    public void testMaxValue() {
        final double maxValue = Double.MAX_VALUE / 1024;
        final double[] sFiniteExponent = new double[] {0,1,-1,-2,0.5,-0.5};
        final double[] sValueExceptionExponent = new double[]{2};
        testSingle( maxValue, 0, sFiniteExponent, sValueExceptionExponent);
        testSingle(-maxValue, 0, sFiniteExponent, sValueExceptionExponent);
        testSingle( maxValue, maxValue, sFiniteExponent, sValueExceptionExponent);
        testSingle(-maxValue, maxValue, sFiniteExponent, sValueExceptionExponent);
        testSingle( maxValue, Double.MIN_NORMAL, sFiniteExponent, sValueExceptionExponent);
        testSingle(-maxValue, Double.MIN_NORMAL, sFiniteExponent, sValueExceptionExponent);
    }

    @Test
    public void testInitLSB() {
        try {
            op = new IntvDbl(1, Double.NaN);
            assertEquals(1, op.value(), TOLERANCE);
            assertEquals(1, op.uncertainty() * (1L << 52), TOLERANCE);

            op = new IntvDbl(Double.MIN_NORMAL, Double.NaN);
            assertEquals(1, op.value() / Double.MIN_NORMAL, TOLERANCE);
            assertEquals(1, op.uncertainty() / Double.MIN_VALUE, TOLERANCE);

            op = new IntvDbl(Double.MIN_NORMAL / 2, Double.NaN);
            assertEquals(1, op.value() / Double.MIN_NORMAL * 2, TOLERANCE);
            assertEquals(1, op.uncertainty() / Double.MIN_VALUE, TOLERANCE);

        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    private void testInitFail(double value) {
        try {
            op = new IntvDbl(value);
            fail();
        } catch (ValueException e) {
        }

        try {
            op = new IntvDbl(value, 0);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }

        try {
            op = new IntvDbl(0, value);
            fail();
        } catch (ValueException e) {
            fail(e.getMessage());
        } catch (UncertaintyException e) {
        }
    }

    @Test
    public void testInitFail() {
        testInitFail(Double.POSITIVE_INFINITY);
        testInitFail(Double.NEGATIVE_INFINITY);

        try {
            op = new IntvDbl(Double.NaN);
        } catch (ValueException e) {
        }

        final double lsbMax = Double.MAX_VALUE * 1E-16;

        try {
            op = new IntvDbl(Double.MAX_VALUE, 0);
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            op = new IntvDbl(Double.MAX_VALUE, lsbMax);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }

        try {
            op = new IntvDbl(0, Double.MAX_VALUE);
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            op = new IntvDbl(-lsbMax, Double.MAX_VALUE);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }
    }

    @Test
    public void testShiftFail() {
        try {
            op = new IntvDbl(Double.MAX_VALUE * 0.999);
            op.shift(1);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            op = new IntvDbl(Double.MAX_VALUE / 4, Double.MAX_VALUE / 2);
            op.shift(2);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }
    }

    @Test
    public void testScaleFail() {
        try {
            op = new IntvDbl(Double.MAX_VALUE * 0.999);
            op.multiply(1.002);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            op = new IntvDbl(Double.MAX_VALUE / 4, Double.MAX_VALUE / 2);
            op.multiply(2.001);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }
    }

    private void testAdd( double value1, double range1, double value2, double range2, boolean isAdd) {
        try {
            op = new IntvDbl(value1, range1);
            res = op.clone();
            op2 = new IntvDbl(isAdd? value2 : -value2, range2);
            res.add(op2);
            final double value = isAdd? value1 + value2 : value1 - value2;
            if (value == 0) {
                assertEquals(0, res.value(), TOLERANCE);
            } else {
                assertEquals(String.format("%s + %s: value %e != %e", op.toString(), op2.toString(), value, res.value()), 
                    1, res.value() / value, TOLERANCE);
            }
            if (Double.isNaN(range1)) {
                range1 = IReal.getLSB(value1);
            }
            if (Double.isNaN(range2)) {
                range2 = IReal.getLSB(value2);
            }
            final double range = range1 + range2;
            assertEquals(String.format("%s + %s: range %e != %e", op.toString(), op2.toString(), range, res.uncertainty()), 
                         range, res.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }
    }
    private void testAdd(final double[] sValue, final double[] sRange) {
        for (double value1: sValue) {
            for (double value2: sValue) {
                for (double range1: sRange) {
                    for (double range2: sRange) {
                        testAdd( value1, range1, value2, range2, true);
                        testAdd( value1, range1, value2, range2, false);
                    }
                }
            }
        }
    }

    @Test
    public void testAdd() {
        testAdd(new double[]{0, Double.MIN_NORMAL, 1, -1, 2, -2}, 
                new double[]{0, 0.1, 10, 1E-3, 1E+3, 1E-6, 1E+6, 1E-9, 1E+9, 1E-12, 1E+12, 1E-15, 1E+15, Double.NaN, Double.MAX_VALUE / 2});
        testAdd(new double[]{0, Double.MIN_VALUE, -Double.MIN_VALUE, Double.MIN_NORMAL, -Double.MIN_NORMAL}, 
                new double[]{0, Double.MIN_VALUE, Double.MIN_NORMAL});
        testAdd(new double[]{0, Double.MIN_VALUE / 2, -Double.MIN_VALUE / 2}, 
                new double[]{0, Double.MIN_VALUE / 2});

        try {
            op = new IntvDbl(Double.MAX_VALUE / 2, Double.MAX_VALUE / 2);
            res = op.clone();
            op2 = new IntvDbl(Double.MAX_VALUE, 0);
            res.add(op2);
            fail();
        } catch (ValueException e) {
        } catch (TypeException | UncertaintyException e) {
            fail(e.getMessage());
        }

        try {
            op = new IntvDbl(0, Double.MAX_VALUE / 2);
            res = op.clone();
            res.add(new IntvDbl(0, Double.MAX_VALUE));
            fail();
        } catch (UncertaintyException e) {
        } catch (TypeException | ValueException e) {
            fail(e.getMessage());
        }

        try {
            op = new IntvDbl();
            res = op.clone();
            res.add(null);
            fail();
        } catch (TypeException e) {
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }

    }


    private void testMultiply(double d1, double d2, boolean isMultiply) {
        final double d = d1 * d2;
        final double tolerance = (isMultiply? TOLERANCE : TOLERANCE_DIV) * Math.abs(d);
        try {
            if ((d + tolerance) < (res.value() - res.uncertainty())) {
                fail(String.format("(%.3e in %s) %s (%.3e in %s): value %.3e < %s by %.3e", 
                    d1, op.toString(), isMultiply? "*" : "/",  d2, 
                    isMultiply? op2.toString() : op2.power(-1).toString(), 
                    d, res.toString(), d - (res.value() - res.uncertainty())));
            }
            if ((res.value() + res.uncertainty()) < (d - tolerance)) {
                fail(String.format("(%.3e in %s) %s (%.3e in %s): value %.3e > %s by %.3e", 
                    d1, op.toString(), isMultiply? "*" : "/",  d2, 
                    isMultiply? op2.toString() : op2.power(-1).toString(), 
                    d, res.toString(), d - (res.value() + res.uncertainty())));
            }
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testMultiply(boolean isMultiply, int samples) {
        final double div = op.uncertainty() / samples;
        final double div2 = op2.uncertainty() / samples;
        if (div > 0) {
            if (((op.value() - op.uncertainty()) <= 0) && 
                (0 <= (op.value() - op.uncertainty()))) {
                    assertTrue((res.value() - res.uncertainty()) <= 0);
                    assertTrue(0 <= (res.value() + res.uncertainty()));
            }
            for (int i = -samples; i <= samples; ++i) {
                final double d1 = op.value() + i*div;
                if (div2 > 0) {
                    for (int j = -samples; j <= samples; ++j) {
                        final double d2 = op2.value() + j*div2;
                        testMultiply(d1, d2, isMultiply);
                    }
                } else {
                    testMultiply(d1, op2.value(), isMultiply);
                }
            }
        } else if (div2 > 0) {
            for (int j = -samples; j <= samples; ++j) {
                final double d2 = op2.value() + j*div2;
                testMultiply(op.value(), d2, isMultiply);
            }
        }
    }
    private void multiply( double value1, double range1, double value2, double range2, boolean isMultiply) throws ValueException, UncertaintyException, TypeException {
        op = new IntvDbl(value1, range1);
        res = op.clone();
        op2 = isMultiply? new IntvDbl(value2, range2) : new IntvDbl(value2, range2).power(-1);
        res.multiply(op2);
    }
    private void testMultiply( double value1, double range1, double value2, double range2, boolean isMultiply) {
        try {
            multiply( value1, range1, value2, range2, isMultiply);
            testMultiply(value1, isMultiply? value2 : 1.0 / value2, isMultiply);
            testMultiply(isMultiply, 8);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testMultiply(final double[] sValue, final double[] sRange) {
        for (double value1: sValue) {
            for (double value2: sValue) {
                for (double range1: sRange) {
                    for (double range2: sRange) {
                        testMultiply( value1, range1, value2, range2, true);
                        if (range2 < Math.abs(value2)) {
                            testMultiply( value1, range1, value2, range2, false);
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testMultiply() {
        final double maxValue = Math.sqrt(Double.MAX_VALUE);
        final double minValue = Math.sqrt(Double.MIN_VALUE);
        final double minNormal = Math.sqrt(Double.MIN_NORMAL);

        testMultiply(new double[]{0, 1, -1, 2, -2}, 
                new double[]{0, 0.1, 10, 1E-3, 1E+3, 1E-6, 1E+6, 1E-9, 1E+9, 1E-12, 1E+12, 1E-15, 1E+15, Double.NaN, maxValue});
        testMultiply(new double[]{0, minValue, -minValue, minNormal, -minNormal}, 
                new double[]{0, minValue, minNormal});

        try {
            op = new IntvDbl(2E-17, 1);
            res = op.clone();
            op2 = new IntvDbl(2E-17, 1);
            res.multiply(op2);
            assertEquals(0, res.value(), 0);
            assertEquals(1, res.uncertainty(), 0);
        } catch (ValueException e) {
        } catch (TypeException | UncertaintyException e) {
            fail(e.getMessage());
        }

        try {
            op = new IntvDbl(maxValue, maxValue);
            res = op.clone();
            op2 = new IntvDbl(maxValue, 0);
            res.multiply(op2);
            fail();
        } catch (UncertaintyException e) {
        } catch (TypeException | ValueException e) {
            fail(e.getMessage());
        }

        try {
            op = new IntvDbl(0, maxValue);
            res = op.clone();
            op2 = new IntvDbl(0, Double.MAX_VALUE);
            res.multiply(op2);
            fail();
        } catch (UncertaintyException e) {
        } catch (TypeException | ValueException e) {
            fail(e.getMessage());
        }

        try {
            op = new IntvDbl();
            res = op.clone();
            res.multiply(null);
            fail();
        } catch (TypeException e) {
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
}

