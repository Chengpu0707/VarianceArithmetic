/*
 * Compare IntvDbl and varDbl on 
 *  *) stablity of mutiplication
 *  *) exagration of subtraction of large mutiplications
 */
package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestDblCompare {
    final double TOLERANCE = 2E-16;

    private void testMultiple(final double value1, final double uncertainty1,
                              final double value2, final double uncertainty2,
                              final double varValue, final double varUncertainty,
                              final double intvValue, final double intvUncertainty,
                              final double intvValueThreshold, final double intvRangeThreshold) {
        try {
            final IntvDbl op1 = new IntvDbl(value1, uncertainty1);
            final IntvDbl op2 = new IntvDbl(value2, uncertainty2);
            final IntvDbl res = op1.multiply(op2);
            assertEquals(intvValue, res.value(), intvValueThreshold);
            assertEquals(intvUncertainty, res.uncertainty(), intvRangeThreshold);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }

        try {
            final VarDbl op1 = new VarDbl(value1, uncertainty1*uncertainty1);
            final VarDbl op2 = new VarDbl(value2, uncertainty2*uncertainty2);
            final VarDbl res = op1.multiply(op2);
            assertEquals(varValue, res.value(), TOLERANCE);
            assertEquals(varUncertainty, res.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    /*
     * variance arithmetic is stable
     * interval arithmetic is not stable
     */
    @Test
    public void testMultiple() {
        testMultiple(1, 0.1, 2, 0.1, 
            2, 0.22383029285599393, 2.01, 0.3, 5E-16, 3E-16);
        testMultiple(1, 0.1, -2, 0.1, 
            -2, 0.22383029285599393, -2.01, 0.3, 
            5E-16, 3E-16);
        testMultiple(-1, 0.1, 2, 0.1, 
            -2, 0.22383029285599393, -2.01, 0.3, 
            5E-16, 3E-16);
        testMultiple(-1, 0.1, -2, 0.1, 
            2, 0.22383029285599393, 2.01, 0.3, 
            5E-16, 3E-16);

            
        testMultiple(1, 0.2, 2, 0.1, 
            2, 0.4127953488110059, 2.02, 0.5, 
            5E-16, 3E-16);
        testMultiple(1, 0.2, -2, 0.1, 
            -2, 0.4127953488110059, -2.02, 0.5, 
            5E-16, 3E-16);
        testMultiple(-1, 0.2, 2, 0.1, 
            -2, 0.4127953488110059, -2.02, 0.5, 
            5E-16, 3E-16);
        testMultiple(-1, 0.2, -2, 0.1, 
            2, 0.4127953488110059, 2.02, 0.5, 
            5E-16, 3E-16);

            
        testMultiple(1, 0.2, 0, 0.1, 
            0, 0.10198039027185571, 0, 0.12, 
            5E-16, 3E-16);
        testMultiple(-1, 0.2, 0, 0.1, 
            0, 0.10198039027185571, 0, 0.12, 
            5E-16, 3E-16);

    }


    @Test
    public void testLargeDiff() {
        assertEquals(13316075197586562L, 64919121L*205117922L);
        assertEquals(13316075197586561L, 159018721L*83739041L);
        assertEquals(1L, 64919121L*205117922L - 159018721L*83739041L);

        testMultiple(64919121, Double.NaN, 205117922, Double.NaN, 
                    6658037598793281L*2, 2.8284271247461903, 6658037598793281L*2, 4, 
                    5E-16, 3E-16); 
        testMultiple(-159018721L, Double.NaN, 83739041, Double.NaN, 
                    -6658037598793280L*2, 3.4641016151377544, -6658037598793280L*2, 5, 
                    5E-16, 3E-16);

        try {
            final VarDbl op1 = new VarDbl( 6658037598793281L*2, 2.8284271247461903*2.8284271247461903);
            final VarDbl op2 = new VarDbl(-6658037598793280L*2, 3.4641016151377544*3.4641016151377544);
            final VarDbl res = op1.add(op2);
            assertEquals(2, res.value(), TOLERANCE);
            assertEquals(4.472135954999579393, res.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException | TypeException e) {
            fail();
        }
        
        try {
            final IntvDbl op1 = new IntvDbl( 6658037598793281L*2, 4);
            final IntvDbl op2 = new IntvDbl(-6658037598793280L*2, 5);
            final IntvDbl res = (IntvDbl) op1.add(op2);
            assertEquals(2, res.value(), TOLERANCE);
            assertEquals(9, res.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException | TypeException e) {
            fail();
        }

    }
    
}
