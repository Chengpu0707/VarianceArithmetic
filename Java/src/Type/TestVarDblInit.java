/*
 * Test initiation of a VarDbl value:
 *  *) When the varince is 0
 */
package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

public class TestVarDblInit extends VarDbl {
	public TestVarDblInit() {
        super();
    }

    private static final double VARIANCE_TOLERANCE = 3E-16;   // 53-bit

    /*
     * special tolerance:
     *  -1: value() == 0
     *  -2: uncertainty = 0
     */
    private void test(final double value, final double dev, final double tolerance) {
        try {
            init(value, dev*dev, false);

            if ((tolerance == -1) ||  (value == 0)) {
                assertTrue(0 == value()); 
            } else {
                assertEquals(1, value() / value, tolerance);
            }

            if ((tolerance == -2) || (dev == 0)) {
                assertTrue(0 == uncertainty()); 
            } else {
                assertEquals(1, uncertainty() / dev, VARIANCE_TOLERANCE);
            }

            if ((tolerance != -1) && (tolerance != -2) &&(0 != value) && (0 < dev)) {
                assertEquals(1, Math.sqrt(precSq()) / (dev / value), tolerance);
            }
        } catch (ValueException e) {
            e.printStackTrace();
            fail();
        } catch (UncertaintyException e) {
            e.printStackTrace();
            fail();
        }
    }
    private void test(final double value, final double dev) {
        test(value, dev, VARIANCE_TOLERANCE);
    }

    @Test
    public void test0() {
        test(0, 0);
        assertEquals("0", toString());
    }

    @Test
    public void test1() {
        test(1, 0);
        assertEquals("1", toString());
    }

    @Test
    public void testNeg1() {
        test(-1, 0);
        assertEquals("-1", toString());
    }

    @Test
    public void testPrec1024() {
        test(1024, 1);
        assertEquals(-26, exp());
        assertEquals(1L << 36, val());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());
        assertEquals("1.024e+03~1", toString());

        test(1, 1.0 / 1024);
        assertEquals(-36, exp());
        assertEquals(1L << 36, val());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());
        assertEquals("1~9.8e-04", toString());
    }

    @Test
    public void testDevVanish() {
        test(1L << 40, 1);
        assertEquals(-12, exp());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, val());
        assertEquals(1 << 24, var());
        assertEquals("1.100e+12~1", toString());

        test(1L << 50, 1);
        assertEquals(-2, exp());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, val());
        assertEquals(1 << 4, var());
        assertEquals("1.126e+15~1", toString());

        test(1L << 52, 1);
        assertEquals(0, exp());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, val());
        assertEquals(1, var());
        assertEquals("4.504e+15~1", toString());

        test(1L << 53, 1, -2);
        assertEquals(1, exp());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, val());
        assertEquals(0, var());
        assertEquals("9.007e+15", toString());

        test(Double.MAX_VALUE, Math.sqrt(Double.MAX_VALUE), -2);
        assertEquals(Dbl.DOUBLE_EXP_MAX, exp());
        assertEquals(Dbl.DOUBLE_VAL_MAX, val());
        assertEquals(0, var());
        assertEquals("1.798e+308", toString());
    }

    @Test
    public void testValVanish() {
        test(1, 1 << 20);
        assertEquals(-6, exp());
        assertEquals(1L << 6, val());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());
        assertEquals("1~1.0e+06", toString());

        test(1, 1 << 26);
        assertEquals(0, exp());
        assertEquals(1, val());
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());
        assertEquals("1~6.7e+07", toString());
        
        test(1, 1 << 27, -1);
        assertEquals(1, exp());
        assertEquals(0, val());     // rounding off
        assertEquals(Dbl.DOUBLE_VAL_EXTRA, var());
        assertEquals("0~1.3e+08", toString());

        test(Double.MIN_NORMAL / 2, Math.sqrt(Double.MIN_VALUE), -1);
        assertEquals(Dbl.DOUBLE_EXP_MIN / 2, exp());
        assertEquals(0, val());
        assertEquals(1, var());
        assertEquals("0~2.2e-162", toString());
    }

    @Test
    public void testExpMin() {
        test(Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE));
        assertEquals(Dbl.DOUBLE_EXP_MIN / 2, exp());
        assertEquals(1L << (Dbl.DOUBLE_EXP_SHIFT / 2), val());
        assertEquals(1, var());
        assertEquals("1.492e-154~2.2e-162", toString());
    }

    @Test
    public void testExpMax() {
        test(Double.MAX_VALUE, Math.sqrt(Double.MAX_VALUE), -2);
        assertEquals(Dbl.DOUBLE_EXP_MAX, exp());
        assertEquals(Dbl.DOUBLE_VAL_MAX, val());
        assertEquals(0, var());
        assertEquals("1.798e+308", toString());

        test(Math.sqrt(Double.MAX_VALUE), Math.sqrt(Double.MAX_VALUE));
        assertEquals((Dbl.DOUBLE_EXP_MAX + 1) / 2, exp());
        assertEquals("1.341e+154~1.3e+154", toString());        
    }


    @Test
    public void testExpOverMax() {
        try {
            init(Math.sqrt(Double.MAX_VALUE), Double.POSITIVE_INFINITY, false);
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            assertEquals("0", toString());
        }
    }

    @Test
    public void testApproxVal() {
        test(1.1, 1E-1, 2E-6);
        assertEquals("1.100e+00~1.0e-01", toString());
        test(1.1, 1E-2, 1E-7);
        assertEquals("1.100e+00~1.0e-02", toString());
        test(1.1, 1E-3, 1E-8);
        assertEquals("1.100e+00~1.0e-03", toString());
        test(1.1, 1E-4, 2E-9);
        assertEquals("1.100e+00~1.0e-04", toString());
        test(1.1, 1E-5, 2E-10);
        assertEquals("1.100e+00~1.0e-05", toString());
        test(1.1, 1E-6, 2E-10);
        assertEquals("1.100e+00~1.0e-06", toString());
    }

    @Test
    public void testInitLSB() {
        VarDbl var;
        try {
            var = new VarDbl(1, Double.NaN);
            assertEquals(1, var.value(), VARIANCE_TOLERANCE);
            assertEquals(1, var.uncertainty() * (1L << 52), VARIANCE_TOLERANCE);

            var = new VarDbl(Double.MIN_NORMAL, Double.NaN);
            assertEquals(1, var.value() / Double.MIN_NORMAL, VARIANCE_TOLERANCE);
            assertEquals(1, var.uncertainty() / Double.MIN_VALUE, VARIANCE_TOLERANCE);

        } catch (ValueException | UncertaintyException e) {
        }

    }

    @Test
    public void testCopy() {
        final VarDbl other = this;
        assertEquals(exp(), other.exp());
        assertEquals(neg(), other.neg());
        assertEquals(val(), other.val());
        assertEquals(var(), other.var());
    }
}
