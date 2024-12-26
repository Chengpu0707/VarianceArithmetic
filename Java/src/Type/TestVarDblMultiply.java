package Type;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class TestVarDblMultiply {

    private enum ETest {
        TEST_ALL, TEST_ZERO_VALUE, TEST_ZERO_DEV,
    };
    
    private void test(final VarDbl varDbl, final double value, final double dev, final double tolerance, ETest eTest) 
            throws InitException {
        if ((value == 0) || (eTest == ETest.TEST_ZERO_VALUE)) {
            assertEquals(0, varDbl.value(), tolerance);
        } else {
            assertEquals(1, varDbl.value() / value, tolerance);
        }

        if ((dev == 0) || (eTest == ETest.TEST_ZERO_DEV)) {
            assertEquals(0, varDbl.uncertainty(), Double.MIN_NORMAL);
        } else {
            assertEquals(dev, varDbl.uncertainty(), Math.ulp(dev));
        }
    }

    private VarDbl test(final double value1, final double dev1, ETest eTest1, 
                        final double value2, final double dev2, ETest eTest2,
                        final double value, final double dev, ETest eTest,
                        final double tolerance) throws InitException {
        final VarDbl op1 = new VarDbl(value1, dev1);
        final VarDbl op2 = new VarDbl(value2, dev2);
        final VarDbl prod1 = op1.multiply(op2);
        test(prod1, value, dev, tolerance, eTest);
        final VarDbl prod2 = op2.multiply(op1);
        test(prod2, value, dev, tolerance, eTest);
        return prod2;
    }
    private VarDbl test(final double value1, final double dev1, 
                        final double value2, final double dev2,
                        final double value, final double dev,
                        final double tolerance) throws InitException {
        return test(value1, dev1, ETest.TEST_ALL, value2, dev2, ETest.TEST_ALL, value, dev, ETest.TEST_ALL, tolerance);
    }

    @Test
    public void testMultiplyPrecise() throws InitException {
        test(1, 0.001, 0, 0, 0, 0, 0);
        test(1, 0.001, 1, 0, 1, 0.001, 0);
        test(1, 0.001, -1, 0, -1, 0.001, 0);
    }

    @Test
    public void testMultiplyZero() throws InitException {
        test(1, 1.0/1024, 0, 0, 0, 0, 0);
        test(1, 1.0/1024, 0, 1.0/1024, 0, Math.sqrt(1.0/(1L << 20) + 1.0/(1L << 40)), 0);
        test(1, 1.0/1024, 0, 2.0/1024, 0, Math.sqrt(4.0/(1L << 20) + 4.0/(1L << 40)), 0);
    }

    @Test
    public void testMultiplyOne() throws InitException {
        test(1, 1.0/1024, 1, 0, 1, 1.0/1024, 0);
        test(1, 1.0/1024, 1, 1.0/1024, 1, Math.sqrt(2.0/(1L << 20) + 1.0/(1L << 40)), 0);
        test(1, 1.0/1024, 1, 2.0/1024, 1, Math.sqrt(5.0/(1L << 20) + 4.0/(1L << 40)), 0);
    }

    @Test
    public void testMultiply() throws InitException {
        test(2, 1.0/1024, 0.5, 0.5/1024, 1, Math.sqrt(1.25/(1L << 20) + 0.25/(1L << 40)), 0);
        test(2, 1.0/1024, 0.5, 1.0/1024, 1, Math.sqrt(4.25/(1L << 20) + 1.0/(1L << 40)), 0);
        test(2, 1.0/1024, 0.5, 2.0/1024, 1, Math.sqrt(16.25/(1L << 20) + 4.0/(1L << 40)), 0);
    }

    @Test
    public void testMultiplyApprox() throws InitException {
        test(2, 1E-3, 0.5, 0.5E-3, 1, Math.sqrt(1.25E-6 + 0.25E-12), 0);
        test(2, 1E-3, 0.5, 1E-3, 1, Math.sqrt(4.25E-6 + 1E-12), 0);
        test(2, 1E-3, 0.5, 2E-3, 1, Math.sqrt(16.25E-6 + 4E-12), 0);
    }

    @Test
    public void testMultiplyMin() throws InitException {
        test(1, 1E-3, Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), 
             Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE + Double.MIN_NORMAL * 1E-6), 0);
        test(Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), Math.sqrt(Double.MIN_NORMAL), Math.sqrt(Double.MIN_VALUE), 
             Double.MIN_NORMAL, 0, 0);
    }

    @Test
    public void testMultiplyMax() throws InitException {
        VarDbl op1 = new VarDbl(1, 1E-6);
        final VarDbl op2 = new VarDbl(Math.sqrt(Double.MAX_VALUE), Math.sqrt(Double.MAX_VALUE / 2));
        final VarDbl prod = op1.multiply(op2);
        assertEquals(Math.sqrt(Double.MAX_VALUE), prod.value(), Math.ulp(prod.value()));
        final double variance = Double.MAX_VALUE * 1E-6 + Double.MAX_VALUE / 2 + 1E-6 * Double.MAX_VALUE / 2;
        assertEquals(1, prod.variance() / variance, 3E-6);
    }

    @Test
    public void testLargeDiff() throws InitException {
        assertEquals(13316075197586562L, 64919121L*205117922L);
        assertEquals(13316075197586561L, 159018721L*83739041L);
        assertEquals(1L, 64919121L*205117922L - 159018721L*83739041L);

        final VarDbl op1 = new VarDbl(64919121L), op2 = new VarDbl(205117922L);
        final VarDbl op3 = new VarDbl(-159018721L), op4 = new VarDbl(83739041L);
        final VarDbl re12 = op1.multiply(op2);
        assertEquals(13316075197586562.0, re12.value(), 0);
        assertEquals(0, re12.uncertainty(), 0);
        final VarDbl re34 = op3.multiply(op4);
        assertEquals( -13316075197586560.0, re34.value(), 0);
        assertEquals(0.5, re34.uncertainty(), 0);
        final VarDbl re = re12.add(re34);
        assertEquals(2, re.value(), Math.ulp(re.value()));
        assertEquals(VarDbl.ulp(re12.value()), VarDbl.ulp(re34.value()), Math.ulp(re.uncertainty()));
        assertEquals(0.5, re.uncertainty(), 0);
    }
}
