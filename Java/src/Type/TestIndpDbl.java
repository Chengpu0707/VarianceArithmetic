package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import Stats.Stat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;

import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestIndpDbl {
    private static final double TOLERANCE = 2E-16;   // 53-bit
    private static final double TOLERANCE_DIV = 4E-16;   // 53-bit

    IndpDbl op;
    IndpDbl op2;
    IndpDbl res;
    
    @Test 
    public void testSimple() {
    }
    @Test
    public void testHash() {
        final Map<IndpDbl, Integer> sMap = new HashMap<IndpDbl, Integer>();
       try {
            final IndpDbl one = new IndpDbl(1, 0);
            final IndpDbl two = new IndpDbl(1, 2);
            final IndpDbl three = new IndpDbl(0, 1);
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

    private void testToString(String str, double value, double dev) {
        try {
            op = new IndpDbl(value, dev*dev);
            assertEquals(str, op.toString());
        } catch (ValueException | UncertaintyException e) {
            fail( e.getMessage());
        }
    }

    @Test
    public void testToString() {
        testToString("0", 0, 0);
        testToString("-1&100", -1, 100);
        testToString("-1&1.0e+03", -1, 1000);    
        testToString("1.000e-03&1.0e-03", 0.001, 0.001);
        testToString("1.000e-03&2.2e-19", 0.001, Double.NaN);
        testToString("-1.341e+154&1.5e+138", -Math.sqrt(Double.MAX_VALUE), Double.NaN);
        testToString("2.225e-308", Double.MIN_NORMAL, Double.MIN_VALUE);
        testToString("2.225e-308&1.3e+154", Double.MIN_NORMAL, Math.sqrt(Double.MAX_VALUE));        
    }


    private void testClone() {
        res = op.clone();
        assertTrue(op.equals(res));
        assertTrue(res.equals(op));
        assertTrue(op != res);
        assertFalse(op == res);

        IndpDbl clone = res.clone();
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
                Stat stat = new Stat();
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
    private void testPower(double value, double dev, double exponent, boolean raiseValueException) {
        try {
            op = new IndpDbl(value, dev*dev);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
        testPower(exponent, raiseValueException);
    }

    private void testSingle(double value, double dev, 
                            final double[] sFiniteExponent, final double[] sValueExceptionExponent) {
        try {
            op = new IndpDbl(value, dev*dev);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
        assertEquals(value, op.value(), TOLERANCE);
        assertEquals(dev, op.uncertainty(), TOLERANCE);
        testClone();
        testNegate();
        testShift();
        testOffset();
        testScale();
        testPower(sFiniteExponent, sValueExceptionExponent);
    }

    
}
