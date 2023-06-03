package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestIndpDbl {
    private static final double VARIANCE_TOLERANCE = 3E-16;   // 53-bit
    
    private enum ETest {
        TEST_ALL, TEST_ZERO_VALUE, TEST_ZERO_DEV,
    };
    
    private void test(final IndpDbl IndpDbl, final double value, final double dev, final double tolerance, ETest eTest) 
            throws ValueException, UncertaintyException {
        if ((value == 0) || (eTest == ETest.TEST_ZERO_VALUE)) {
            assertEquals(0, IndpDbl.value(), tolerance);
        } else {
            assertEquals(1, IndpDbl.value() / value, tolerance);
        }

        if ((dev == 0) || (eTest == ETest.TEST_ZERO_DEV)) {
            assertEquals(0, IndpDbl.uncertainty(), VARIANCE_TOLERANCE);
        } else {
            assertEquals(1, IndpDbl.uncertainty() / dev, VARIANCE_TOLERANCE);
        }
    }

    private IndpDbl test(final double value1, final double dev1, ETest eTest1, 
                        final double value2, final double dev2, ETest eTest2,
                        final double value, final double dev, ETest eTest,
                        final double tolerance) {
        try {
            final IndpDbl op1 = new IndpDbl(value1, dev1*dev1);
            final IndpDbl op2 = new IndpDbl(value2, dev2*dev2);
            final IndpDbl prod1 = (IndpDbl) op1.multiply(op2);
            test(prod1, value, dev, tolerance, eTest);
            final IndpDbl prod2 = (IndpDbl) op2.multiply(op1);
            test(prod2, value, dev, tolerance, eTest);
            return prod2;
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
            return null;
        }
    }
    private IndpDbl test(final double value1, final double dev1, 
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
        test(1, 1.0/1024, 0, 1.0/1024, 0, Math.sqrt(1.0/(1L << 20)), 1E-16);
        test(1, 1.0/1024, 0, 2.0/1024, 0, Math.sqrt(4.0/(1L << 20)), 1E-16);
    }

    @Test
    public void testMultiplyOne() {
        test(1, 1.0/1024, 1, 0, 1, 1.0/1024, 1E-16);
        test(1, 1.0/1024, 1, 1.0/1024, 1, Math.sqrt(2.0/(1L << 20)), 1E-16);
        test(1, 1.0/1024, 1, 2.0/1024, 1, Math.sqrt(5.0/(1L << 20)), 1E-16);
    }

    @Test
    public void testMultiply() {
        test(2, 1.0/1024, 0.5, 0.5/1024, 1, Math.sqrt(1.25/(1L << 20)), 3E-16);
        test(2, 1.0/1024, 0.5, 1.0/1024, 1, Math.sqrt(4.25/(1L << 20)), 3E-16);
        test(2, 1.0/1024, 0.5, 2.0/1024, 1, Math.sqrt(16.25/(1L << 20)), 3E-16);
    }

    @Test
    public void testMultiplyApprox() {
        test(2, 1E-3, 0.5, 0.5E-3, 1, Math.sqrt(1.25E-6), 3E-16);
        test(2, 1E-3, 0.5, 1E-3, 1, Math.sqrt(4.25E-6), 3E-16);
        test(2, 1E-3, 0.5, 2E-3, 1, Math.sqrt(16.25E-6), 3E-16);
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
        IndpDbl op1;
        try {
            op1 = new IndpDbl(1, 1E-6);
            final IndpDbl op2 = new IndpDbl(Math.sqrt(Double.MAX_VALUE), Double.MAX_VALUE / 2);
            final IndpDbl prod = (IndpDbl) op1.multiply(op2);
            assertEquals(1, prod.value() / Math.sqrt(Double.MAX_VALUE), 3E-16);
            final double variance = Double.MAX_VALUE * 1E-6 + Double.MAX_VALUE / 2;
            assertEquals(1, prod.variance() / variance, 3E-16);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testLargeDiff() {
        assertEquals(13316075197586562L, 64919121L*205117922L);
        assertEquals(13316075197586561L, 159018721L*83739041L);
        assertEquals(1L, 64919121L*205117922L - 159018721L*83739041L);

        try {
            final IndpDbl op1 = new IndpDbl(64919121L), op2 = new IndpDbl(205117922L);
            final IndpDbl op3 = new IndpDbl(-159018721L), op4 = new IndpDbl(83739041L);
            final IndpDbl re12 = (IndpDbl) op1.multiply(op2);
            assertEquals(6658037598793281L, re12.val());
            assertEquals(1, re12.exp());
            assertEquals(2, re12.var());
            final IndpDbl re34 = (IndpDbl) op3.multiply(op4);
            assertEquals(6658037598793280L, re34.val());
            assertEquals(1, re34.exp());
            assertEquals(3, re34.var());
            final IndpDbl re = (IndpDbl) re12.add(re34);
            assertEquals(2L << 24, re.val());
            assertEquals(-24, re.exp());
            assertEquals(5L << 50, re.var());
            assertEquals(2, re.value(), 1E-16);
            assertEquals(2 * Math.sqrt(5), re.uncertainty(), 1E-16);
        } catch (ValueException | UncertaintyException | TypeException e) {
            fail();
        }
    }
    
    private void testPower( double value1, double dev1, double exponent, 
                            double value, double dev, double tolerance ) {
        try {
            final IndpDbl op1 = new IndpDbl(value1, dev1*dev1);
            final IndpDbl intv = (IndpDbl) op1.power(exponent);
            assertEquals(value, intv.value(), tolerance);
            assertEquals(dev, intv.uncertainty(), tolerance);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    private void testPower( double value1, double dev1, double exponent) {

        try {
            final IndpDbl op1 = new IndpDbl(value1, dev1);
            op1.power(exponent);
            fail();
        } catch (ValueException | UncertaintyException e) {
        }
    }
        
    @Test
    public void testPower() {
        testPower(0, 0.1, 0, 1, 0, VARIANCE_TOLERANCE);
        testPower(1, 0.1, 0, 1, 0, VARIANCE_TOLERANCE);
        testPower(2, 0.1, 0, 1, 0, VARIANCE_TOLERANCE);
        testPower(-1, 0.1, 0, 1, 0, VARIANCE_TOLERANCE);
        testPower(1, 1, 0, 1, 0, VARIANCE_TOLERANCE);
        testPower(1, Double.NaN, 0, 1, 0, VARIANCE_TOLERANCE);

        testPower(0, 0.1, 1, 0, 0.1, VARIANCE_TOLERANCE);
        testPower(1, 0.1, 1, 1, 0.1, VARIANCE_TOLERANCE);
        testPower(2, 0.1, 1, 2, 0.1, VARIANCE_TOLERANCE);
        testPower(-1, 0.1, 1, -1, 0.1, VARIANCE_TOLERANCE);
        testPower(1, 1, 1, 1, 1, VARIANCE_TOLERANCE);
        testPower(1, Double.NaN, 0, 1, 1.0/(1L << Dbl.DOUBLE_EXP_SHIFT), VARIANCE_TOLERANCE);

        testPower(0, 0.1, 2, 0, 0, VARIANCE_TOLERANCE);
        testPower(1, 0.1, 2, 1, 0.2, VARIANCE_TOLERANCE);
        testPower(2, 0.1, 2, 4, 0.8, VARIANCE_TOLERANCE);
        testPower(-1, 0.1, 2, 1, 0.2, VARIANCE_TOLERANCE);
        testPower(1, 1, 2, 1, 2, VARIANCE_TOLERANCE);
        testPower(1, Double.NaN, 2, 1, 2.0/(1L << Dbl.DOUBLE_EXP_SHIFT), VARIANCE_TOLERANCE);

        testPower(0, 0.1, 0.5, 0, 0, VARIANCE_TOLERANCE);
        testPower(1, 0.1, 0.5, 1, 0.05, VARIANCE_TOLERANCE);
        testPower(2, 0.1, 0.5, 1.4142135623842478, 0.05 * Math.sqrt(2), VARIANCE_TOLERANCE);
        testPower(-1, 0.1, 0.5);
        testPower(1, 1, 0.5, 1, 0.5, VARIANCE_TOLERANCE);
        testPower(1, Double.NaN, 0, 1, 0.5/(1L << Dbl.DOUBLE_EXP_SHIFT), VARIANCE_TOLERANCE);

        testPower(0, 0.1, -1);
        testPower(1, 0.1, -1, 1, 0.1, VARIANCE_TOLERANCE);
        testPower(2, 0.1, -1, 0.5, 0.05, VARIANCE_TOLERANCE);
        testPower(-1, 0.1, -1, -1, 0.1, VARIANCE_TOLERANCE);
        testPower(1, 1, -1, 1, 1, VARIANCE_TOLERANCE );
        testPower(1, Double.NaN, 0, 1, 1.0/(1L << Dbl.DOUBLE_EXP_SHIFT), VARIANCE_TOLERANCE);
    }
    
}
