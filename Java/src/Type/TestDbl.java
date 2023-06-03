package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.ValueException;

public class TestDbl {
	private static final double TOLERANCE = 1E-16;

    private void test(double d, long val, long exp) {
        final Dbl dbl = new Dbl(d);
        try {
            assertEquals(val, dbl.val);
            assertEquals(exp, dbl.exp);
            assertEquals(d, dbl.toDouble(), TOLERANCE);
        } catch (ValueException e) {
            e.printStackTrace();
            fail();
        }
    }

    private void test(int exp, boolean neg, long val, double d) {
        final Dbl dbl = new Dbl(exp, neg, val);
        try {
            assertEquals(d, dbl.toDouble(), TOLERANCE);
       } catch (ValueException e) {
            e.printStackTrace();
            fail();
        }
    }

    private void testException(long val, int exp) {
        final Dbl dbl = new Dbl(exp, true, val);
        try {
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
    public void testMax() {
        test(Double.MAX_VALUE, Dbl.DOUBLE_VAL_MAX, Dbl.DOUBLE_EXP_MAX);
    }

    @Test
    public void testMin() {
        test(Double.MIN_VALUE, 1, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinHalf() {
        test(Double.MIN_VALUE / 2, 0, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinDoble() {
        test(Double.MIN_VALUE * 2, 2, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormal() {
        test(Double.MIN_NORMAL, Dbl.DOUBLE_VAL_EXTRA, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormalHalf() {
        test(Double.MIN_NORMAL / 2, Dbl.DOUBLE_VAL_EXTRA / 2, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormalThreeQuater() {
        test(Double.MIN_NORMAL * 3 / 4, Dbl.DOUBLE_VAL_EXTRA * 3 / 4, Dbl.DOUBLE_EXP_MIN);
    }

    @Test
    public void testMinNormalSqrt() {
        test(Math.sqrt(Double.MIN_NORMAL), Dbl.DOUBLE_VAL_EXTRA, Dbl.DOUBLE_EXP_MIN  / 2 - 26);
    }

    @Test
    public void testPositiveInfinity() {
        assertFalse(Double.isNaN(Double.POSITIVE_INFINITY));
        assertTrue(Double.isInfinite(Double.POSITIVE_INFINITY));
        test(Double.POSITIVE_INFINITY, Dbl.DOUBLE_VAL_EXTRA, Dbl.DOUBLE_EXP_MAX + 1);
    }

    @Test
    public void testNegativeInfinity() {
        assertFalse(Double.isNaN(Double.NEGATIVE_INFINITY));
        assertTrue(Double.isInfinite(Double.NEGATIVE_INFINITY));
        test(Double.NEGATIVE_INFINITY, Dbl.DOUBLE_VAL_EXTRA, Dbl.DOUBLE_EXP_MAX + 1);
    }

    @Test
    public void testNaN() {
        assertTrue(Double.isNaN(Double.NaN));
        assertFalse(Double.isInfinite(Double.NaN));
        test(Double.NaN, Dbl.DOUBLE_VAL_NAN, Dbl.DOUBLE_EXP_MAX + 1);
    }

    @Test
    public void testAssign() {
        test(- Dbl.DOUBLE_EXP_SHIFT + 4, true, Dbl.DOUBLE_VAL_EXTRA, -16);
    }

    @Test
    public void testNormalization() {
        test(0, true, 1, -1);
    }

    @Test
    public void testMoreThanMaxVal() {
        testException(Dbl.DOUBLE_VAL_MAX + 1, Dbl.DOUBLE_EXP_MAX);
    }

    @Test
    public void testMoreThanMaxExp() {
        testException(Dbl.DOUBLE_VAL_MAX, Dbl.DOUBLE_EXP_MAX + 1);
    }

    @Test
    public void testLessThanMinVal() {
        test(Dbl.DOUBLE_EXP_MIN, false, 1, 0);
    }

    @Test
    public void testLessThanMinExp() {
        test(Dbl.DOUBLE_EXP_MIN - 1, false, 1, 0);
    }

   @Test
    public void testNegativeValMin() {
        testException(-1, 0);
    }

    @Test
    public void testDividZero() {
        assertEquals(Double.NEGATIVE_INFINITY, -1.0 / 0, TOLERANCE);

        try {
            int i = 1 / 0;
            fail();
        } catch( Exception e) {
        }
    }

    @Test
    public void testSqrtNegative() {
        assertEquals(Double.NaN, Math.sqrt( -1 ), TOLERANCE);
    }
    
    @Test
    public void testSqrtInfinity() {
        assertEquals(Double.POSITIVE_INFINITY, Math.sqrt( Double.POSITIVE_INFINITY ), TOLERANCE);
    }

    @Test
    public void testRoundingError() {
        assertEquals(1, Math.sqrt(2) * Math.sqrt(2) / 2, 3E-16);
        assertEquals(1, Math.sqrt(2) / Math.sqrt(0.5) / 2, 3E-16);
        assertEquals(1, Math.sqrt(0.5) * Math.sqrt(0.5) / 0.5, 3E-16);

        assertEquals(1, Math.sqrt(5) * Math.sqrt(5) / 5, 3E-16);

        assertEquals(1, 1.0/3 + 2.0/3, 1E-16);
        assertEquals(1, 1.0/11 + 10.0/11, 1E-16);
    }
}
