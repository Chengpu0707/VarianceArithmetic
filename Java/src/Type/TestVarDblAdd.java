package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestVarDblAdd {
    private static final double VARIANCE_TOLERANCE = 5E-16;   // 53-bit
    
    private void test(final VarDbl varDbl, final double value, final double dev, final double tolerance) 
            throws ValueException, UncertaintyException {
        if (value == 0) {
            assertEquals(0, varDbl.value(), tolerance);
        } else {
            assertEquals(1, varDbl.value() / value, tolerance);
        }

        if (dev == 0) {
            assertEquals(0, varDbl.uncertainty(), VARIANCE_TOLERANCE);
        } else {
            assertEquals(1, varDbl.uncertainty() / dev, VARIANCE_TOLERANCE);
        }
    }

    private VarDbl test(final double value1, final double dev1, 
                        final double value2, final double dev2,
                        final double value, final double dev,
                        final double tolerance) {
        try {
            final VarDbl op1 = new VarDbl(value1, dev1*dev1);
            final VarDbl op2 = new VarDbl(value2, dev2*dev2);
            final VarDbl sum1 = op1.clone();
            sum1.add(op2);
            test(sum1, value, dev, tolerance);
            final VarDbl sum2 = op2.clone();
            sum2.add(op1);
            test(sum2, value, dev, tolerance);
            return sum2;
        } catch (ValueException | UncertaintyException e) {
            fail();
            return null;
        }
    }

    @Test
    public void testAddPrecise() {
        test(1, 0.001, 0, 0, 1, 0.001, 1E-16);
        test(1, 0.001, 1, 0, 2, 0.001, 1E-16);
        test(1, 0.001, -1, 0, 0, 0.001, 1E-16);
    }

    @Test
    public void testAddZero() {
        test(1, 1.0/1024, 0, 0, 1, 1.0/1024, 1E-16);
        test(1, 1.0/1024, 0, 1.0/1024, 1, Math.sqrt(2)/1024, 1E-16);
        test(1, 1.0/1024, 0, 2.0/1024, 1, Math.sqrt(5)/1024, 1E-16);
    }

    @Test
    public void testAdd() {
        test(1.0/3, 1.0/1024, 2.0/3, 1.0/ (1L << 52), 1, 1.0/1024, 1E-16);
        test(1.0/3, 1.0/1024, 2.0/3, 1.0/1024, 1, Math.sqrt(2)/1024, 1E-16);
        test(1.0/3, 1.0/1024, 2.0/3, 2.0/1024, 1, Math.sqrt(5)/1024, 1E-16);
    }

    @Test
    public void testAddApprox() {
        test(1.0/3, 0.001, 2.0/3, 1E-16, 1, 0.001, 1E-16);
        test(1.0/3, 0.001, 2.0/3, 0.001, 1, Math.sqrt(2) * 0.001, 1E-16);
        test(1.0/3, 0.001, 2.0/3, 0.002, 1, Math.sqrt(5) * 0.001, 1E-16);
        test(1.0/3, 0.001, 2.0/3, 0.01, 1, Math.sqrt(1E-4 + 1E-6), 1E-16);
    }

    @Test
    public void testAddRoundingError() {
        try {
            final VarDbl op1 = new VarDbl(1.0/3, 1.0/(1L << 20));
            final VarDbl op2 = new VarDbl(2.0/3, 4.0/(1L << 20));
            assertEquals(22906492245L, op1.val());
            assertEquals(22906492245L, op2.val());
            assertEquals(-36, op1.exp());
            assertEquals(-35, op2.exp());
            assertEquals(1.0 * op1.val() /(1L << 36), op1.value(), 1E-16);
            assertEquals(-4.8506E-12, op1.value() - 1.0/3, 1E-16);
            assertEquals(-9.7012E-12, op2.value() - 2.0/3, 1E-16);
            // The rounding error accumulates when it is always rounded off in addition
            assertEquals(-1.4552E-11, op1.value() + op2.value() - 1.0, 1E-16);
            assertEquals(-2.91038E-11, 1.0 * (op1.val() /2 + op2.val()) /(1L << 35) - 1.0, 1E-16);
            // Whe the rounding error is accounted for
            assertEquals(true, op1.rndv());
            assertEquals(true, op2.rndv());
            assertEquals(0, 1.0 * (1 + op1.val() /2 + op2.val()) /(1L << 35) - 1.0, 1E-16);
            test(1.0/3, 0.001, 2.0/3, 0.002, 1, Math.sqrt(5) * 0.001, 1E-16);

            // When the rounding errors do not cancel each other, conventional double calc gives better result
            final Dbl res = new Dbl(op1.value() - op2.value());
            assertEquals(4.8506E-12, res.toDouble() + 1.0/3, 1E-16);
            final long val = 1 + op1.val() /2 - op2.val();
            assertEquals(1.94025E-11, 1.0 * val /(1L << 35) + 1.0/3, 1E-16);
            final VarDbl dif = (VarDbl) op1.add(op2.negate());
            assertEquals(-35, dif.exp());
            assertEquals(Math.abs(val), dif.val());
            assertEquals(1.94025E-11, dif.value() + 1.0/3, 1E-16);
//            assertEquals(0, dif.add(new VarDbl(1.0/3)).value(), 1E-16);

        } catch (ValueException | UncertaintyException e) {
           fail();;
        }
    }

    @Test
    public void testSbtraction() {
        test(1.0/3, 1.0/1024, -2.0/3, 0, -1.0/3, 1.0/1024, 3E-10);
        test(1.0/3, 0, -2.0/3, 1.0/1024, -1.0/3, 1.0/1024, 3E-10);
        test(1.0/3, 1.0/1024, -2.0/3, 1.0/1024, -1.0/3, Math.sqrt(2)/1024, 3E-10);
        test(1.0/3, 2.0/1024, -2.0/3, 1.0/1024, -1.0/3, Math.sqrt(5)/1024, 3E-10);
        test(1.0/3, 1.0/1024, -2.0/3, 2.0/1024, -1.0/3, Math.sqrt(5)/1024, 3E-10);

        test(1, 1.0/1024, -1, 2.0/1024, 0, Math.sqrt(5)/1024, 1E-16);

        final double dev = Math.sqrt(0.001*0.001 + 0.002*0.002);
        test(1.001, 0.001, -1, 0.002, 0.001, dev, 2E-7);
    }

    @Test
    public void testValVanish() {
        VarDbl sum;
        sum = test(1.0/3, 1 << 20, 2.0/3, 1 << 20, 1, Math.sqrt(2) * (1 << 20), 1E-16);
        assertEquals(32, sum.val());
        sum = test(1.0/3, 1 << 20, 2.0/3, 2 << 20, 1, Math.sqrt(5) * (1 << 20), 1E-16);
        assertEquals(32, sum.val());
        sum = test(1.0/3, 2 << 20, 2.0/3, 2 << 20, 1, Math.sqrt(2) * (2 << 20), 1E-16);
        assertEquals(16, sum.val());
    }

    @Test
    public void testAddMin() {
        test(1, 1.0/1024, Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), 
             1, 1.0/1024, 1E-16);
        test(Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), 
             Math.sqrt(Double.MIN_NORMAL) * 2, Math.sqrt(Double.MIN_VALUE * 2), 1E-16);
    }

    @Test
    public void testAddMax() {
        test(1, 1.0/1024, Math.sqrt(Double.MAX_VALUE), Math.sqrt(Double.MAX_VALUE / 2), 
             Math.sqrt(Double.MAX_VALUE), Math.sqrt(Double.MAX_VALUE / 2), 3E-16);
        test(Math.sqrt(Double.MAX_VALUE) / 2, Math.sqrt(Double.MAX_VALUE / 4), 
             Math.sqrt(Double.MAX_VALUE) / 4, Math.sqrt(Double.MAX_VALUE / 4), 
             Math.sqrt(Double.MAX_VALUE) * 3 / 4, Math.sqrt(Double.MAX_VALUE / 2), 3E-16);
    }

}
