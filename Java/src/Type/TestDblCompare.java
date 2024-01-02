/*
 * Compare IntvDbl and varDbl on 
 *  *) stablity of mutiplication
 *  *) exagration of subtraction of large mutiplications
 */
package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestDblCompare {

    private void testMultiple(final double value1, final double uncertainty1,
                              final double value2, final double uncertainty2,
                              final double varValue, final double varUncertainty,
                              final double intvValue, final double intvUncertainty) {
        try {
            final IntvDbl op1 = new IntvDbl(value1, uncertainty1);
            final IntvDbl op2 = new IntvDbl(value2, uncertainty2);
            final IntvDbl res = op1.multiply(op2);
            assertEquals(intvValue, res.value(), 5E-16);
            assertEquals(intvUncertainty, res.uncertainty(), 3E-16);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }

        try {
            final VarDbl op1 = new VarDbl(value1, uncertainty1);
            final VarDbl op2 = new VarDbl(value2, uncertainty2);
            final VarDbl res = op1.multiply(op2);
            assertEquals(varValue, res.value(), 2E-16);
            assertEquals(varUncertainty, res.uncertainty(), 2E-16);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }
    private void testMultipleAll(final double value1, final double uncertainty1,
                              final double value2, final double uncertainty2,
                              final double varValue, final double varUncertainty,
                              final double intvValue, final double intvUncertainty) {
        testMultiple(value1, uncertainty1, value2, uncertainty2,
                     varValue, varUncertainty, intvValue, intvUncertainty);
        testMultiple(-value1, uncertainty1, value2, uncertainty2,
                     -varValue, varUncertainty, -intvValue, intvUncertainty);
        testMultiple(value1, uncertainty1, -value2, uncertainty2,
                     -varValue, varUncertainty, -intvValue, intvUncertainty);
        testMultiple(-value1, uncertainty1, -value2, uncertainty2,
                     varValue, varUncertainty, intvValue, intvUncertainty);
    }

    /*
     * Both variance arithmetic and interval arithmetic are stable
     */
    @Test
    public void testMultiple() {
        testMultipleAll(1, 0.1, 2, 0.1, 
            2, 0.22383029285599393, 
            2.01, 0.3000000000000007);
        testMultipleAll(1, 0.2, 2, 0.1, 
            2, 0.4127953488110059, 
            2.02, 0.5000000000000004);
        testMultipleAll(1, 0.2, 0, 0.1, 
            0, 0.10198039027185571, 
            0, 0.12);
        testMultipleAll(0.1, 0.2, 0.2, 0.1, 
            0.02, 0.045825756949558406, 
            0.03, 0.06);
    }


    @Test
    public void testLargeDiff() {
        assertEquals(13316075197586562L, 64919121L * 205117922L);
        assertEquals(13316075197586561L, 159018721L * 83739041L);
        assertEquals(1L, 64919121L * 205117922L - 159018721L * 83739041L);

        testMultiple(64919121, 0, 205117922, 0, 
                     13316075197586562L, 1.1547005383792517, 
                     13316075197586562L, 2); 
        testMultiple(-159018721L, 0, 83739041, 0, 
                    -13316075197586560L, 1.1547005383792517, 
                    -13316075197586560L, 2);

        try {
            final VarDbl op1 = new VarDbl( 13316075197586562L, 1.1547005383792517);
            final VarDbl op2 = new VarDbl(-13316075197586560L, 1.1547005383792517);
            final VarDbl res = op1.add(op2);
            assertEquals(2, res.value(), res.uncertainty());
            assertEquals(1.6329931618554523, res.uncertainty(), 
                         Dbl.getLSB(res.uncertainty()));
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
        
        try {
            final IntvDbl op1 = new IntvDbl( 13316075197586562L, 2);
            final IntvDbl op2 = new IntvDbl(-13316075197586560L, 2);
            final IntvDbl res = (IntvDbl) op1.add(op2);
            assertEquals(2, res.value(), res.uncertainty());
            assertEquals(4, res.uncertainty(), Dbl.getLSB(res.uncertainty()));
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

    }
    
}
