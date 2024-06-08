package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;

import java.util.HashMap;
import java.util.Map;


public class TestVarDbl {
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

        } catch (InitException e) {
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
        } catch (InitException e) {
            fail(e.getMessage());
        }
    }

    private void testToString(String str, double value, double dev) {
        try {
            op = new VarDbl(value, dev);
            assertEquals(str, op.toString());
        } catch (InitException e) {
            fail(e.getMessage());
        }
     }
 
    private void testToString(String str, double value) {
        try {
            op = new VarDbl(value);
            assertEquals(str, op.toString());
        } catch (InitException e) {
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
        testToString("1.341e+154", Math.sqrt(Double.MAX_VALUE));
        testToString("2.225e-308", Double.MIN_NORMAL, Double.MIN_VALUE);
        testToString("2.225e-308", Double.MIN_NORMAL);        
    }
 
    @Test 
    public void testInitLong() {
        op = new VarDbl(VarDbl.DOUBLE_MAX_SIGNIFICAND - 1);
        assertEquals(VarDbl.DOUBLE_MAX_SIGNIFICAND - 1, op.value(), 0);
        assertEquals(0, op.uncertainty(), 0);

        op = new VarDbl(VarDbl.DOUBLE_MAX_SIGNIFICAND);
        assertEquals(VarDbl.DOUBLE_MAX_SIGNIFICAND, op.value(), 0);
        assertEquals(0, op.uncertainty(), 0);

        op = new VarDbl(VarDbl.DOUBLE_MAX_SIGNIFICAND + 1);
        assertEquals(VarDbl.DOUBLE_MAX_SIGNIFICAND + 1, op.value(), 0);
        assertEquals(0.5, op.uncertainty(), 0);

        op = new VarDbl(VarDbl.DOUBLE_MAX_SIGNIFICAND + 2);
        assertEquals(VarDbl.DOUBLE_MAX_SIGNIFICAND + 2, op.value(), 0);
        assertEquals(0, op.uncertainty(), 0);
    }

    @Test
    public void testInitMax() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2);
        } catch( InitException e) {
            fail(e.getMessage());
        }
    }
 
    private void testClone() {
        res = op.clone();
        assertFalse(op.equals(res));
        assertFalse(res.equals(op));
        assertTrue(op != res);
        assertFalse(op == res);
        try {
            assertEquals(0, op.compareTo(res));
        } catch (InitException e) {
            fail(e.getMessage());
        }
        try {
            assertEquals(0, res.compareTo(op));
        } catch (InitException e) {
            fail(e.getMessage());
        }

        VarDbl clone = res.clone();
        assertFalse(op.equals(clone));
        assertFalse(clone.equals(op));
        try {
            assertEquals(0, op.compareTo(clone));
        } catch (InitException e) {
            fail(e.getMessage());
        }
        try {
            assertEquals(0, clone.compareTo(op));
        } catch (InitException e) {
            fail(e.getMessage());
        }
    }
    
    private void testNegate() {
        res = op.clone();
        res.negate();
        assertEquals(-op.value(), res.value(), 0);
        final double uncertainty = res.uncertainty();
        assertEquals(op.uncertainty(), uncertainty, 0);
    }
    
    private void testOffset( double offset ) {
        try {
            op2 = new VarDbl(offset);
            res = op.add(op2);
            assertEquals(op.value() + offset, res.value(), 0);
            assertEquals(Math.sqrt(op.uncertainty()*op.uncertainty() + op2.uncertainty()*op2.uncertainty()), 
                         res.uncertainty(), Math.ulp(res.uncertainty()));
        } catch (InitException e) {
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
            op2 = new VarDbl(fold);
            res = op.multiply(op2);
            assertEquals(op.value() * fold, res.value(), res.ulp());
            if (op.uncertainty() > 0) {
                double exp = op.uncertainty() * Math.abs(fold);
                double lsb = op2.uncertainty() * Math.abs(op.value());
                if (exp < lsb)
                    exp = lsb;
                final double uncertainty = res.uncertainty();
                assertEquals(exp, uncertainty, Math.ulp(uncertainty));
            }
        } catch (InitException e) {
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
            assertEquals(value, op.value(), op.ulp());
            assertEquals(Math.sqrt(dev*dev), op.uncertainty(), Math.ulp(op.uncertainty()));
        } catch (InitException e) {
            fail(e.getMessage());
        }
        testClone();
        testNegate();
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
    public void testScaleOverflow() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2, Math.sqrt(Double.MAX_VALUE/2));
            res.multiply(new VarDbl(4));
            fail();
        } catch (InitException e) {
        }
    }

    @Test
    public void testOffsetOverflow() {
        try {
            res = new VarDbl(Double.MAX_VALUE/2, Math.sqrt(Double.MAX_VALUE/2));
            res.add(new VarDbl(Double.MAX_VALUE));
            fail();
        } catch (InitException e) {
         }
    }

}
