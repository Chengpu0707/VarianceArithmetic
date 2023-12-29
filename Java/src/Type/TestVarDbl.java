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
    VarDbl op;
    VarDbl op2;
    VarDbl res;
    
    @Test 
    public void testSimple() {
        try {
            double value = 0;
            double dev =  Double.MIN_NORMAL;
            op = new VarDbl(value, dev*dev);
            final double tolerance = Math.pow(2, op.exp());
            assertEquals(value, op.value(), tolerance);
            assertEquals(Math.sqrt(dev*dev), op.uncertainty(), tolerance);
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
            op = new VarDbl(value, dev*dev);
            assertEquals(str, op.toString());
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
     }
 
    @Test
    public void testToString() {
        testToString("0", 0, 0);
        testToString("-1~100", -1, 100);
        testToString("-1~1.0e+03", -1, 1000);    
        testToString("1.000e-03~1.0e-03", 0.001, 0.001);
        testToString("1.000e-03~2.2e-19", 0.001,Double.NaN);
        testToString("1.341e+154~1.5e+138", Math.sqrt(Double.MAX_VALUE), Double.NaN);
        testToString("2.225e-308", Double.MIN_NORMAL, Double.MIN_VALUE);
        testToString("2.225e-308", Double.MIN_NORMAL, Double.NaN);        
    }
 
    @Test
    public void testInitMax() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2, Double.NaN);
        } catch( UncertaintyException ex) {
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
        try {
            res = op.clone();
            res.negate();
            final double tolerance = Math.pow(2, res.exp());
            assertEquals(-op.value(), res.value(), tolerance);
            assertEquals(op.uncertainty(), res.uncertainty(), tolerance);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    
    private void testShift( int bits ) {
        try {
            res = op.clone();
            res.shift(bits);
            final double tolerance = Math.pow(2, res.exp());
            if (bits >= 0) {
                assertEquals(op.value() * (1L << bits), res.value(), tolerance);
                assertEquals(op.uncertainty() * (1L << bits), res.uncertainty(), tolerance);
            } else {
                assertEquals(op.value() / (1L << -bits), res.value(), tolerance);
                assertEquals(op.uncertainty() / (1L << -bits), res.uncertainty(), tolerance);
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
            final double tolerance = Math.pow(2, res.exp());
            assertEquals(op.value() + offset, res.value(), tolerance);
            assertEquals(op.uncertainty(), res.uncertainty(), tolerance);
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
            res.multiply(fold);
            final double tolerance = Math.pow(2, res.exp());
            assertEquals(op.value() * fold, res.value(), tolerance);
            assertEquals(op.uncertainty() * Math.abs(fold), res.uncertainty(), tolerance);
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
            op = new VarDbl(value, dev*dev);
            final double tolerance = Math.pow(2, op.exp());
            assertEquals(value, op.value(), tolerance);
            assertEquals(Math.sqrt(dev*dev), op.uncertainty(), tolerance);
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
            res.multiply(4);
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
            res.add(Double.MAX_VALUE);
            fail();
        } catch (ValueException e) {
        } catch (UncertaintyException e) {
            fail(e.getMessage());
        }
    }

}
