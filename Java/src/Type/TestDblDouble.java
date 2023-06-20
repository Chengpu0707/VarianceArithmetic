package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;


import Type.IReal.ValueException;

/*
 * Test conversion between double and Dbl.
 */
public class TestDblDouble {
	private static final double TOLERANCE = 1E-16;

    private void test(double d, long val, long exp) {
        Dbl dbl;
        try {
            dbl = new Dbl(d);
            assertEquals(val, dbl.val);
            assertEquals(exp, dbl.exp);
            assertEquals(d, dbl.toDouble(), TOLERANCE);
        } catch (ValueException e) {
            fail();
        }
    }

    private void test(int exp, boolean neg, long val, double d) {
        try {
            final Dbl dbl = new Dbl(exp, neg, val, false);
            assertEquals(d, dbl.toDouble(), TOLERANCE);
       } catch (ValueException e) {
            e.printStackTrace();
            fail();
        }
    }

    private void testException(double d) {
        try {
            new Dbl(d);
            fail();
        } catch (ValueException e) {
        }
    }

    private void testException(int exp, long val) {
        try {
            Dbl dbl = new Dbl(exp, false, val, false);
            dbl.toDouble();
            fail();
        } catch (ValueException e) {
        }
    }

    @Test
    public void test0() {
        test(0, 0, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void test1() {
        test(1, Dbl.DOUBLE_VAL_EXTRA, - Dbl.DOUBLE_EXP_SHIFT);
    }

    @Test
    public void testNeg1() {
        test(-1, Dbl.DOUBLE_VAL_EXTRA, - Dbl.DOUBLE_EXP_SHIFT);;
    }
    
    @Test
    public void test16() {
        test(16, Dbl.DOUBLE_VAL_EXTRA, - Dbl.DOUBLE_EXP_SHIFT + 4);
    }

    @Test
    public void testAssign() {
        test(- Dbl.DOUBLE_EXP_SHIFT + 4, true, Dbl.DOUBLE_VAL_EXTRA, -16);
    }

    @Test
    public void testMax() {
        test(Double.MAX_VALUE, Dbl.DOUBLE_VAL_MAX, Dbl.DOUBLE_EXP_MAX);
    }

    @Test
    public void testMin() {
        test(Double.MIN_VALUE, 1, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinHalved() {
        test(Double.MIN_VALUE / 2, 0, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinDoubled() {
        test(Double.MIN_VALUE * 2, 2, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormal() {
        test(Double.MIN_NORMAL, Dbl.DOUBLE_VAL_EXTRA, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormalHalfed() {
        test(Double.MIN_NORMAL / 2, Dbl.DOUBLE_VAL_EXTRA / 2, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormalSqrt() {
        test(Math.sqrt(Double.MIN_NORMAL), Dbl.DOUBLE_VAL_EXTRA, Dbl.DOUBLE_EXP_MIN  / 2 - Dbl.DOUBLE_EXP_SHIFT / 2);
    }

    @Test
    public void testNotFiniteInit() {
        testException(Double.POSITIVE_INFINITY);
        testException(Double.NEGATIVE_INFINITY);
        testException(Double.NaN);
    }

    @Test
    public void testNegativeVal() {
        testException(0, -1);
        testException(0, Long.MAX_VALUE + 1);
    }

    @Test
    public void testNormalization() {
        test(0, true, 1, -1);
        test(0, false, 1, 1);
        test(0, true, Dbl.DOUBLE_VAL_EXTRA << 1, -(1L << 53));
        test(0, false, Dbl.DOUBLE_VAL_EXTRA << 1, (1L << 53));
    }

    @Test
    public void testMoreThanMaxVal() {
        testException(Dbl.DOUBLE_EXP_MAX, Dbl.DOUBLE_VAL_MAX + 1L);
    }

    @Test
    public void testMoreThanMaxExp() {
        testException(Dbl.DOUBLE_EXP_MAX + 1, Dbl.DOUBLE_VAL_EXTRA);
    }

    @Test
    public void testLessThanMinVal() {
        test(Dbl.DOUBLE_EXP_MIN, false, 1, 0);
    }

    @Test
    public void testLessThanMinExp() {
        test(Dbl.DOUBLE_EXP_MIN - 1, false, 1, 0);
    }
}
