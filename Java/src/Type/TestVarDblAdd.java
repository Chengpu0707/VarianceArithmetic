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
            final VarDbl op1 = new VarDbl(value1, dev1);
            final VarDbl op2 = new VarDbl(value2, dev2);
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
    public void testAddLSB() {
        test(1.0, Math.ulp(1.0), 1.0, Math.ulp(1.0), 2.0, 0.0, 1E-16);
        test(1.0/3, Math.ulp(1.0/3), 2.0/3, Math.ulp(2.0/3), 1.0, 0.0, 1E-16);
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
