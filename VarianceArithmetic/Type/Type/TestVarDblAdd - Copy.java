package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestVarDblAdd {
    private static final double VARIANCE_TOLERANCE = 1E-9;   // 31-bit
    
    private enum ETest {
        TEST_ALL, TEST_ZERO_VALUE, TEST_ZERO_DEV,
    };
    
    private void test(final VarDbl varDbl, final double value, final double dev, final double tolerance, ETest eTest) 
            throws ValueException, UncertaintyException {
        if ((value == 0) || (eTest == ETest.TEST_ZERO_VALUE)) {
            assertEquals(0, varDbl.value(), tolerance);
        } else {
            assertEquals(1, varDbl.value() / value, tolerance);
        }

        if ((dev == 0) || (eTest == ETest.TEST_ZERO_DEV)) {
            assertEquals(0, varDbl.uncertainty(), VARIANCE_TOLERANCE);
        } else {
            assertEquals(1, varDbl.uncertainty() / dev, VARIANCE_TOLERANCE);
        }
    }

    private void test(final double value1, final double dev1, ETest eTest1, 
                      final double value2, final double dev2, ETest eTest2,
                      final double value, final double dev, ETest eTest,
                      final double tolerance) {
        try {
            final VarDbl op1 = new VarDbl(value1, dev1*dev1);
            final VarDbl op2 = new VarDbl(value2, dev2*dev2);
            final VarDbl sum1 = (VarDbl) op1.add(op2);
            test(op1, value1, dev1, tolerance, eTest1);
            test(op2, value2, dev2, tolerance, eTest2);
            test(sum1, value, dev, tolerance, eTest);
            final VarDbl sum2 = (VarDbl) op2.add(op1);
            test(op1, value1, dev1, tolerance, eTest1);
            test(op2, value2, dev2, tolerance, eTest2);
            test(sum2, value, dev, tolerance, eTest);
        } catch (TypeException e) {
            fail();
        } catch (ValueException e) {
            fail();
        } catch (UncertaintyException e) {
            fail();
        }
    }
    private void test(final double value1, final double dev1, 
                       final double value2, final double dev2,
                       final double value, final double dev,
                       final double tolerance) {
        test(value1, dev1, ETest.TEST_ALL, value2, dev2, ETest.TEST_ALL, value, dev, ETest.TEST_ALL, tolerance);
    }

    @Test
    public void testAddZero() {
        test(1, 1.0/1024, 0, 0, 1, 1.0/1024, 1E-16);
        test(1, 1.0/1024, 0, 1.0/1024, 1, Math.sqrt(2)/1024, 1E-16);
        test(1, 1.0/1024, 0, 2.0/1024, 1, Math.sqrt(5)/1024, 1E-16);
    }

    @Test
    public void testAdd() {
        test(1.0/3, 1.0/1024, 2.0/3, 0, 1, 1.0/1024, 3E-8);
        test(1.0/3, 1.0/1024, 2.0/3, 1.0/1024, 1, Math.sqrt(2)/1024, 3E-8);
        test(1.0/3, 1.0/1024, 2.0/3, 2.0/1024, 1, Math.sqrt(5)/1024, 3E-8);
    }

    @Test
    public void testSbtraction() {
        test(1.0/3, 1.0/1024, -2.0/3, 0, -1.0/3, 1.0/1024, 1E-7);
        test(1.0/3, 0, -2.0/3, 1.0/1024, -1.0/3, 1.0/1024, 1E-7);
        test(1.0/3, 1.0/1024, -2.0/3, 1.0/1024, -1.0/3, Math.sqrt(2)/1024, 1E-7);
        test(1.0/3, 2.0/1024, -2.0/3, 1.0/1024, -1.0/3, Math.sqrt(5)/1024, 1E-7);
        test(1.0/3, 1.0/1024, -2.0/3, 2.0/1024, -1.0/3, Math.sqrt(5)/1024, 1E-7);

        test(1, 1.0/1024, -1, 2.0/1024, 0, Math.sqrt(5)/1024, 1E-16);
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
             Math.sqrt(Double.MAX_VALUE), Math.sqrt(Double.MAX_VALUE), 3E-16);
        test(Math.sqrt(Double.MAX_VALUE) / 2, Math.sqrt(Double.MAX_VALUE / 4), 
             Math.sqrt(Double.MAX_VALUE) / 4, Math.sqrt(Double.MAX_VALUE / 4), 
             Math.sqrt(Double.MAX_VALUE) * 3 / 4, Math.sqrt(Double.MAX_VALUE / 2), 3E-16);
    }

}
