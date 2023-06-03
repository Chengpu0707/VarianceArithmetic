package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestIntvDbl {
    IntvDbl intv;
    final double TOLERANCE = 5E-16;
    
    private void testNegate() {
        IntvDbl neg = (IntvDbl) intv.negate();
        try {
            assertEquals(-intv.value(), neg.value(), TOLERANCE);
            assertEquals(intv.uncertainty(), neg.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }
    
    private void testShift() {
        testShift(0);
        testShift(10);
        testShift(-10);
    }
    private void testShift( int bits ) {
        IntvDbl sh;
        try {
            sh = (IntvDbl) intv.shift(bits);
            if (bits >= 0) {
                assertEquals(intv.value() * (1L << bits), sh.value(), TOLERANCE);
                assertEquals(intv.uncertainty() * (1L << bits), sh.uncertainty(), TOLERANCE);
            } else {
                assertEquals(intv.value() / (1L << -bits), sh.value(), TOLERANCE);
                assertEquals(intv.uncertainty() / (1L << -bits), sh.uncertainty(), TOLERANCE);
            }
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    private void testScale() {
        testScale(0);
        testScale(1000);
        testScale(-1000);
        testScale(0.001);
        testScale(-0.001);
    }
    private void testScale( double fold ) {
        IntvDbl sh;
        try {
            sh = (IntvDbl) intv.scale(fold);
            assertEquals(intv.value() * fold, sh.value(), TOLERANCE);
            assertEquals(intv.uncertainty() * Math.abs(fold), sh.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testInit() {
        try {
            intv = new IntvDbl();
            assertEquals(0, intv.value(), TOLERANCE);
            assertEquals(0, intv.uncertainty(), TOLERANCE);
            testNegate();
            testShift();
            testScale();

            intv = new IntvDbl(1);
            assertEquals(1, intv.value(), TOLERANCE);
            assertEquals(0, intv.uncertainty(), TOLERANCE);
            testNegate();
            testShift();
            testScale();

            intv = new IntvDbl(Double.MAX_VALUE / 2);
            assertEquals(Double.MAX_VALUE / 2, intv.value(), TOLERANCE);
            assertEquals(1, intv.uncertainty() / Double.MAX_VALUE * (1L << 54), TOLERANCE);
            testNegate();

            intv = new IntvDbl(0, 1);
            assertEquals(0, intv.value(), TOLERANCE);
            assertEquals(1, intv.uncertainty(), TOLERANCE);
            testNegate();
            testShift();
            testScale();

            intv = new IntvDbl(1, -1);
            assertEquals(1, intv.value(), TOLERANCE);
            assertEquals(1, intv.uncertainty(), TOLERANCE);
            testNegate();
            testShift();
            testScale();

            intv = new IntvDbl(Double.MAX_VALUE / 2, Double.MAX_VALUE / 2);
            assertEquals(Double.MAX_VALUE / 2, intv.value(), TOLERANCE);
            assertEquals(Double.MAX_VALUE / 2, intv.uncertainty(), TOLERANCE);
            testNegate();

            intv = new IntvDbl(-Double.MAX_VALUE / 2, -Double.MAX_VALUE / 2);
            assertEquals(-Double.MAX_VALUE / 2, intv.value(), TOLERANCE);
            assertEquals(Double.MAX_VALUE / 2, intv.uncertainty(), TOLERANCE);
            testNegate();

        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testInitFail() {
        try {
            intv = new IntvDbl(Double.NaN);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            intv = new IntvDbl(Double.MAX_VALUE);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }

        try {
            intv = new IntvDbl(Double.NEGATIVE_INFINITY);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            intv = new IntvDbl(1, Double.POSITIVE_INFINITY);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }

        try {
            intv = new IntvDbl(1, Double.NEGATIVE_INFINITY);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }

        try {
            intv = new IntvDbl(Double.MAX_VALUE / 2, Double.MAX_VALUE);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }

        try {
            intv = new IntvDbl(-Double.MAX_VALUE / 2, Double.MAX_VALUE);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }
    }

    @Test
    public void testInitLSB() {
        try {
            intv = new IntvDbl(1, Double.NaN);
            assertEquals(1, intv.value(), TOLERANCE);
            assertEquals(1, intv.uncertainty() * (1L << 52), TOLERANCE);

            intv = new IntvDbl(Double.MIN_NORMAL, Double.NaN);
            assertEquals(1, intv.value() / Double.MIN_NORMAL, TOLERANCE);
            assertEquals(1, intv.uncertainty() / Double.MIN_VALUE, TOLERANCE);

            intv = new IntvDbl(Double.MIN_NORMAL / 2, Double.NaN);
            assertEquals(1, intv.value() / Double.MIN_NORMAL * 2, TOLERANCE);
            assertEquals(1, intv.uncertainty() / Double.MIN_VALUE, TOLERANCE);

        } catch (ValueException | UncertaintyException e) {
        }

    }


    @Test
    public void testShiftFail() {
        try {
            intv = new IntvDbl(Double.MAX_VALUE * 0.999);
            intv.shift(1);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            intv = new IntvDbl(Double.MAX_VALUE / 4, Double.MAX_VALUE / 2);
            intv.shift(2);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }
    }

    @Test
    public void testScaleFail() {
        try {
            intv = new IntvDbl(Double.MAX_VALUE * 0.999);
            intv.scale(1.002);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail();
        }

        try {
            intv = new IntvDbl(Double.MAX_VALUE / 4, Double.MAX_VALUE / 2);
            intv.scale(2.001);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
        }
    }


    @Test
    public void testAdd() {
        try {
            final IntvDbl op1 = new IntvDbl(1);
            final IntvDbl op2 = new IntvDbl(2);;
            intv = (IntvDbl) op1.add(op2);
            assertEquals(3, intv.value(), TOLERANCE);
            assertEquals(2.0 / (1L << 52), intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(1);
            final IntvDbl op2 = new IntvDbl(-2, 0.01);
            intv = (IntvDbl) op1.add(op2);
            assertEquals(-1, intv.value(), TOLERANCE);
            assertEquals(0.01, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(1, 0.001);
            final IntvDbl op2 = new IntvDbl(1, 0.001);;
            intv = (IntvDbl) op1.add(op2);
            assertEquals(2, intv.value(), TOLERANCE);
            assertEquals(0.002, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(1, 0.001);
            final IntvDbl op2 = new IntvDbl(-1, 0.001);;
            intv = (IntvDbl) op1.add(op2);
            assertEquals(0, intv.value(), TOLERANCE);
            assertEquals(0.002, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(1, 0.001);
            final IntvDbl op2 = new IntvDbl(1, 0.01);;
            intv = (IntvDbl) op1.add(op2);
            assertEquals(2, intv.value(), TOLERANCE);
            assertEquals(0.011, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(Double.MIN_VALUE, Double.MAX_VALUE / 2);
            final IntvDbl op2 = new IntvDbl(-Double.MIN_NORMAL, Double.MAX_VALUE / 2);;
            intv = (IntvDbl) op1.add(op2);
            assertEquals(-Double.MIN_NORMAL + Double.MIN_VALUE, intv.value(), TOLERANCE);
            assertEquals(Double.MAX_VALUE, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(Double.MIN_VALUE, Double.MAX_VALUE / 2);
            final IntvDbl op2 = new IntvDbl(-Double.MIN_NORMAL, Double.MAX_VALUE);;
            intv = (IntvDbl) op1.add(op2);
            fail();
        } catch (UncertaintyException e) {
        } catch (TypeException | ValueException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(Double.MAX_VALUE / 2, Double.MAX_VALUE / 2);
            final IntvDbl op2 = new IntvDbl(Double.MAX_VALUE * 3 / 4, Double.MAX_VALUE / 4);;
            intv = (IntvDbl) op1.add(op2);
            fail();
        } catch (ValueException e) {
        } catch (TypeException | UncertaintyException e) {
            fail();
        }
    }


    @Test
    public void testMultiple() {
        try {
            final IntvDbl op1 = new IntvDbl(1);
            final IntvDbl op2 = new IntvDbl(-2);
            intv = (IntvDbl) op1.multiply(op2);
            assertEquals(-2, intv.value(), TOLERANCE);
            assertEquals(2.0 / (1L << 52), intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final IntvDbl op1 = new IntvDbl(1, 0.1);
            final IntvDbl op2 = new IntvDbl(2, 0.1);
            intv = (IntvDbl) op1.multiply(op2);
            assertEquals(2.01, intv.value(), TOLERANCE*2);
            assertEquals(0.3, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }
        
        try {
            final IntvDbl op1 = new IntvDbl(1, 0.1);
            final IntvDbl op2 = new IntvDbl(0, 0.1);
            intv = (IntvDbl) op1.multiply(op2);
            assertEquals(0, intv.value(), TOLERANCE);
            assertEquals(0.11, intv.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }
    }

    private void testPower( double value1, double range1, double exponent, 
                            double value, double range, double tolerance ) {
        try {
            final IntvDbl op1 = new IntvDbl(value1, range1);
            intv = (IntvDbl) op1.power(exponent);
            assertEquals(value, intv.value(), tolerance);
            assertEquals(range, intv.uncertainty(), tolerance);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    private void testPower( double value1, double range1, double exponent) {

        try {
            final IntvDbl op1 = new IntvDbl(value1, range1);
            op1.power(exponent);
            fail();
        } catch (ValueException | UncertaintyException e) {
        }
    }
        
    @Test
    public void testPower() {
        testPower(0, 0.1, 0, 1, 0, TOLERANCE);
        testPower(1, 0.1, 0, 1, 0, TOLERANCE);
        testPower(2, 0.1, 0, 1, 0, TOLERANCE);
        testPower(-1, 0.1, 0, 1, 0, TOLERANCE);
        testPower(1, 1, 0, 1, 0, TOLERANCE);
        testPower(1, Double.NaN, 0, 1, 0, TOLERANCE);

        testPower(0, 0.1, 1, 0, 0.1, TOLERANCE);
        testPower(1, 0.1, 1, 1, 0.1, TOLERANCE);
        testPower(2, 0.1, 1, 2, 0.1, TOLERANCE);
        testPower(-1, 0.1, 1, -1, 0.1, TOLERANCE);
        testPower(1, 1, 1, 1, 1, TOLERANCE);
        testPower(1, Double.NaN, 0, 1, 1.0/(1L << Dbl.DOUBLE_EXP_SHIFT), TOLERANCE);

        testPower(0, 0.1, 2, 0.005, 0.005, TOLERANCE);
        testPower(1, 0.1, 2, 1.01, 0.2, TOLERANCE);
        testPower(2, 0.1, 2, 4.01, 0.4, TOLERANCE);
        testPower(-1, 0.1, 2, 1.01, 0.2, TOLERANCE);
        testPower(1, 1, 2, 2, 2, TOLERANCE);
        testPower(1, Double.NaN, 0, 1, 2.0/(1L << Dbl.DOUBLE_EXP_SHIFT), TOLERANCE);

        testPower(0, 0.1, 0.5);
        testPower(1, 0.1, 0.5, 0.99874607311033267, 0.05006277505981887, TOLERANCE);
        testPower(2, 0.1, 0.5, 1.41377127491398302, 0.03536639970496084, TOLERANCE);
        testPower(-1, 0.1, 0.5);
        testPower(1, 1, 0.5, 0.7071067811865475, 0.7071067811865475, TOLERANCE);
        testPower(1, Double.NaN, 0.5, 1, 0.5/(1L << Dbl.DOUBLE_EXP_SHIFT), TOLERANCE);

        testPower(0, 0.1, -1);
        testPower(1, 0.1, -1, 1.0101010101010101, 0.10101010101010101, TOLERANCE);
        testPower(2, 0.1, -1, 0.5012531328320802, 0.02506265664160401, TOLERANCE);
        testPower(-1, 0.1, -1, -1.0101010101010101, 0.10101010101010101, TOLERANCE);
        testPower(1, 1, -1);
        testPower(1, Double.NaN, 0, 1, 1.0/(1L << Dbl.DOUBLE_EXP_SHIFT), TOLERANCE);
    }
    
}
