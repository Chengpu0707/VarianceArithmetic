package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;

import java.util.HashMap;
import java.util.Map;

import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestVarDbl {
    static final double tolerance = Double.MIN_NORMAL;
    VarDbl op;
    VarDbl op2;
    VarDbl res;
    
    @Test
    public void testCompareTo() {
        try {
            op = new VarDbl(1.001, 0.001);
            op2 = new VarDbl(1.000, 0.002);
            assertEquals(0, op.compareTo(op2));
            assertEquals(1, op.compareTo(op2, 0.4));
        
            op = new VarDbl(1.002, 0.001);
            assertEquals(1, op.compareTo(op2));
            assertEquals(-1, op2.compareTo(op));

        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }


    @Test
    public void testHash() {
        final Map<VarDbl, Integer> sMap = new HashMap<VarDbl, Integer>();
        try {
            final VarDbl one = new VarDbl(1, 0);
            final VarDbl two = new VarDbl(1, 2);
            final VarDbl three = new VarDbl(0, 1);
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
            op = new VarDbl(value, dev);
            assertEquals(str, op.toString());
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
     }
 
    private void testToString(String str, double value) {
        try {
            op = new VarDbl(value);
            assertEquals(str, op.toString());
        } catch (ValueException e) {
            fail(e.getMessage());
        }
     }
 
    @Test
    public void testToString() {
        testToString("0", 0, 0);
        testToString("-1~100", -1, 100);
        testToString("-1~1.0e+03", -1, 1000);    
        testToString("1.000e-03~1.0e-03", 0.001, 0.001);
        testToString("1.000e-03~1.3e-19", 0.001);
        testToString("1.341e+154~8.6e+137", Math.sqrt(Double.MAX_VALUE));
        testToString("2.225e-308", Double.MIN_NORMAL, Double.MIN_VALUE);
        testToString("2.225e-308", Double.MIN_NORMAL);        
    }
 
    @Test 
    public void testInitLong() {
        op = new VarDbl(1);
        assertEquals(1, op.value(), tolerance);
        assertEquals(0, op.uncertainty(), tolerance);

        op = new VarDbl(Long.MAX_VALUE);
        assertEquals(9.223372036854776E18, op.value(), 1E17);
        assertEquals(2048 * VarDbl.DEVIATION_OF_LSB, op.uncertainty(), tolerance);
    }

    @Test
    public void testInitMax() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2);
        } catch( ValueException e) {
            fail(e.getMessage());
        }
    }
 
     private void testClone() {
        res = op.clone();
        assertTrue(op.equals(res));
        assertTrue(res.equals(op));
        assertTrue(op != res);
        assertFalse(op == res);

        VarDbl clone = res.clone();
        assertTrue(op.equals(clone));
        assertTrue(clone.equals(op));
    }
    
    private void testNegate() {
        res = op.clone();
        res.negate();
        assertEquals(-op.value(), res.value(), 
                     (res.value() == 0)? tolerance : res.getLSB());
        final double uncertainty = res.uncertainty();
        assertEquals(op.uncertainty(), uncertainty, 
                     (uncertainty == 0)? tolerance : Dbl.getLSB(uncertainty));
    }
    
    private void testShift( int bits ) {
        try {
            res = op.clone();
            res.shift(bits);
            if (bits >= 0) {
                assertEquals(op.value() * (1L << bits), res.value(), 
                             (res.value() == 0)? tolerance : res.getLSB());
                final double uncertainty = res.uncertainty();
                assertEquals(op.uncertainty() * (1L << bits), uncertainty, 
                             (uncertainty == 0)? tolerance : Dbl.getLSB(uncertainty));
            } else {
                assertEquals(op.value() / (1L << -bits), res.value(), 
                             (res.value() == 0)? tolerance : res.getLSB());
                final double uncertainty = res.uncertainty();
                assertEquals(op.uncertainty() / (1L << -bits), uncertainty, 
                             (uncertainty == 0)? tolerance : Dbl.getLSB(uncertainty));
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
            res.add(new VarDbl(offset));
            assertEquals(op.value() + offset, res.value(), 
                         (res.value() == 0)? tolerance : res.getLSB());
            if (op.uncertainty() > 0) {
                final double uncertainty = res.uncertainty();
                assertEquals(op.uncertainty(), uncertainty, Dbl.getLSB(uncertainty));
            }
        } catch (ValueException | UncertaintyException e) {
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
            op2 = new VarDbl(fold);
            res.multiply(op2);
            assertEquals(op.value() * fold, res.value(), 
                         (res.value() == 0)? tolerance : res.getLSB());
            if (op.uncertainty() > 0) {
                double exp = op.uncertainty() * Math.abs(fold);
                double lsb = op2.uncertainty() * Math.abs(op.value());
                if (exp < lsb)
                    exp = lsb;
                final double uncertainty = res.uncertainty();
                assertEquals(exp, uncertainty, Dbl.getLSB(uncertainty));
            }
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

    private void testSingle(double value, double dev) {
        try {
            op = new VarDbl(value, dev);
            assertEquals(value, op.value(), 
                         (op.value() == 0)? tolerance : op.getLSB());
            final double uncertainty = op.uncertainty();
            assertEquals(Math.sqrt(dev*dev), uncertainty, 
                         (uncertainty == 0)? tolerance : op.getLSB());
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
        testClone();
        testNegate();
        testShift();
        testOffset();
        testScale();
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
        testSingle( 0, Math.sqrt(Double.MAX_VALUE) / (1L << 10));
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
    public void test10() {
        testSingle(10, 0);
        testSingle(-10, 0);
        testSingle(10, 1);
        testSingle(-10, 1);
        testSingle(10, 10);
        testSingle(-10, 10);
        testSingle(10, 0.1);
        testSingle(-10, 0.1);
        testSingle(10, 1000);
        testSingle(-10, 1000);
        testSingle(10, 0.001);
        testSingle(-10, 0.001);
    }

    @Test
    public void test10th() {
        testSingle(0.1, 0);
        testSingle(-0.1, 0);
        testSingle(0.1, 1);
        testSingle(-0.1, 1);
        testSingle(0.1, 10);
        testSingle(-0.1, 10);
        testSingle(0.1, 0.1);
        testSingle(-0.1, 0.1);
        testSingle(0.1, 1000);
        testSingle(-0.1, 1000);
        testSingle(0.1, 0.001);
        testSingle(-0.1, 0.001);
    }

    @Test
    public void testMinNormal() {
        testSingle( Double.MIN_NORMAL, 0);
        testSingle(-Double.MIN_NORMAL, 0);
        testSingle( Double.MIN_NORMAL, Double.MIN_VALUE);
        testSingle(-Double.MIN_NORMAL, Double.MIN_VALUE);
        testSingle( Double.MIN_NORMAL, Double.MIN_NORMAL);
        testSingle(-Double.MIN_NORMAL, Double.MIN_NORMAL);
    }

    @Test
    public void testMinValue() {
        testSingle( Double.MIN_VALUE, 0);
        testSingle(-Double.MIN_VALUE, 0);
        testSingle( Double.MIN_VALUE, Double.MIN_VALUE);
        testSingle(-Double.MIN_VALUE, Double.MIN_VALUE);
        testSingle( Double.MIN_VALUE, Double.MIN_NORMAL);
        testSingle(-Double.MIN_VALUE, Double.MIN_NORMAL);
    }

    @Test
    public void testMaxValue() {
        final double maxValue = Math.sqrt(Double.MAX_VALUE) / 1024;
        testSingle( maxValue, 0);
        testSingle(-maxValue, 0);
        testSingle( maxValue, maxValue);
        testSingle(-maxValue, maxValue);
        // the half bit dominants
        testSingle( maxValue, Math.sqrt(maxValue));
        testSingle(-maxValue, Math.sqrt(maxValue));
    }

    @Test
    public void testShiftOverflow() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2, Math.sqrt(Double.MAX_VALUE/2));
            res.shift(2);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }

    }

    @Test
    public void testScaleOverflow() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2, Math.sqrt(Double.MAX_VALUE/2));
            res.multiply(new VarDbl(4));
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testOffsetOverflow() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2, Math.sqrt(Double.MAX_VALUE/2));
            res.add(new VarDbl(Double.MAX_VALUE));
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }
    }

}
