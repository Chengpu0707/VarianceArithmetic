package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.ValueException;
import Type.IReal.UncertaintyException;

public class TestVarDblMultiply {
    private static final double VARIANCE_TOLERANCE = 3E-16;   // 53-bit
    
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

    private VarDbl test(final double value1, final double dev1, ETest eTest1, 
                        final double value2, final double dev2, ETest eTest2,
                        final double value, final double dev, ETest eTest,
                        final double tolerance) {
        try {
            final VarDbl op1 = new VarDbl(value1, dev1*dev1);
            final VarDbl op2 = new VarDbl(value2, dev2*dev2);
            final VarDbl prod1 = op1.clone();
            prod1.multiply(op2);
            test(prod1, value, dev, tolerance, eTest);
            final VarDbl prod2 = op2.clone();
            prod2.multiply(op1);
            test(prod2, value, dev, tolerance, eTest);
            return prod2;
        } catch (ValueException | UncertaintyException e) {
            fail();
            return null;
        }
    }
    private VarDbl test(final double value1, final double dev1, 
                        final double value2, final double dev2,
                        final double value, final double dev,
                        final double tolerance) {
        return test(value1, dev1, ETest.TEST_ALL, value2, dev2, ETest.TEST_ALL, value, dev, ETest.TEST_ALL, tolerance);
    }

    @Test
    public void testMultiplyPrecise() {
        test(1, 0.001, 0, 0, 0, 0, 1E-16);
        test(1, 0.001, 1, 0, 1, 0.001, 1E-16);
        test(1, 0.001, -1, 0, -1, 0.001, 1E-16);
    }

    @Test
    public void testMultiplyZero() {
        test(1, 1.0/1024, 0, 0, 0, 0, 1E-16);
        test(1, 1.0/1024, 0, 1.0/1024, 0, Math.sqrt(1.0/(1L << 20) + 1.0/(1L << 40)), 1E-16);
        test(1, 1.0/1024, 0, 2.0/1024, 0, Math.sqrt(4.0/(1L << 20) + 4.0/(1L << 40)), 1E-16);
    }

    @Test
    public void testMultiplyOne() {
        test(1, 1.0/1024, 1, 0, 1, 1.0/1024, 1E-16);
        test(1, 1.0/1024, 1, 1.0/1024, 1, Math.sqrt(2.0/(1L << 20) + 1.0/(1L << 40)), 1E-16);
        test(1, 1.0/1024, 1, 2.0/1024, 1, Math.sqrt(5.0/(1L << 20) + 4.0/(1L << 40)), 1E-16);
    }

    @Test
    public void testMultiply() {
        test(2, 1.0/1024, 0.5, 0.5/1024, 1, Math.sqrt(1.25/(1L << 20) + 0.25/(1L << 40)), 3E-16);
        test(2, 1.0/1024, 0.5, 1.0/1024, 1, Math.sqrt(4.25/(1L << 20) + 1.0/(1L << 40)), 3E-16);
        test(2, 1.0/1024, 0.5, 2.0/1024, 1, Math.sqrt(16.25/(1L << 20) + 4.0/(1L << 40)), 3E-16);
    }

    @Test
    public void testMultiplyApprox() {
        test(2, 1E-3, 0.5, 0.5E-3, 1, Math.sqrt(1.25E-6 + 0.25E-12), 3E-16);
        test(2, 1E-3, 0.5, 1E-3, 1, Math.sqrt(4.25E-6 + 1E-12), 3E-16);
        test(2, 1E-3, 0.5, 2E-3, 1, Math.sqrt(16.25E-6 + 4E-12), 3E-16);
    }

    @Test
    public void testMultiplyMin() {
        test(1, 1E-3, Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), 
             Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE + Double.MIN_NORMAL * 1E-6), 1E-16);
        test(Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), 
             0, 0, 1E-16);
    }

    @Test
    public void testMultiplyMax() {
        VarDbl op1;
        try {
            op1 = new VarDbl(1, 1E-6);
            final VarDbl op2 = new VarDbl(Math.sqrt(Double.MAX_VALUE), Double.MAX_VALUE / 2);
            final VarDbl prod = (VarDbl) op1.multiply(op2);
            assertEquals(1, prod.value() / Math.sqrt(Double.MAX_VALUE), 3E-16);
            final double variance = Double.MAX_VALUE * 1E-6 + Double.MAX_VALUE / 2 + 1E-6 * Double.MAX_VALUE / 2;
            assertEquals(1, prod.variance() / variance, 3E-16);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testLargeDiff() {
        assertEquals(13316075197586562L, 64919121L*205117922L);
        assertEquals(13316075197586561L, 159018721L*83739041L);
        assertEquals(1L, 64919121L*205117922L - 159018721L*83739041L);

        try {
            final VarDbl op1 = new VarDbl(64919121L), op2 = new VarDbl(205117922L);
            final VarDbl op3 = new VarDbl(-159018721L), op4 = new VarDbl(83739041L);
            final VarDbl re12 = op1.clone();
            re12.multiply(op2);
            assertEquals(6658037598793281L, re12.val());
            assertEquals(1, re12.exp());
            assertEquals(2, re12.var());
            final VarDbl re34 = op3.clone();
            re34.multiply(op4);
            assertEquals(6658037598793280L, re34.val());
            assertEquals(1, re34.exp());
            assertEquals(3, re34.var());
            final VarDbl re = re12.clone();
            re.add(re34);
            assertEquals(2L << 24, re.val());
            assertEquals(-24, re.exp());
            assertEquals(5L << 50, re.var());
            assertEquals(2, re.value(), 1E-16);
            assertEquals(2 * Math.sqrt(5), re.uncertainty(), 1E-16);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
        

    }
}
