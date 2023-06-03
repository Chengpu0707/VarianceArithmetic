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
    final double TOLERANCE = 3E-16;

    private void testMultiple(final double value1, final double uncertainty1,
                              final double value2, final double uncertainty2,
                              final double varValue, final double varUncertainty,
                              final double intvValue, final double intvUncertainty) {
        try {
            final IntvDbl op1 = new IntvDbl(value1, uncertainty1);
            final IntvDbl op2 = new IntvDbl(value2, uncertainty2);
            final IntvDbl res = (IntvDbl) op1.multiply(op2);
            assertEquals(intvValue, res.value(), TOLERANCE*2);
            assertEquals(intvUncertainty, res.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }

        try {
            final VarDbl op1 = new VarDbl(value1, uncertainty1*uncertainty1);
            final VarDbl op2 = new VarDbl(value2, uncertainty2*uncertainty2);
            final VarDbl res = (VarDbl) op1.multiply(op2);
            assertEquals(varValue, res.value(), TOLERANCE);
            assertEquals(varUncertainty, res.uncertainty(), TOLERANCE);
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail();
        }
    }

    /*
     * variance arithmetic is stable
     * interval arithmetic is not stable
     */
    @Test
    public void testMultiple() {
        testMultiple(1, 0.1, 2, 0.1, 
            2, 0.22383029285599393, 2.01, 0.3);
        testMultiple(1, 0.1, -2, 0.1, 
            -2, 0.22383029285599393, -2.1000000000000005, 0.21);
        testMultiple(-1, 0.1, 2, 0.1, 
            -2, 0.22383029285599393, -1.9, 0.18999999999999995);
        testMultiple(-1, 0.1, -2, 0.1, 
            2, 0.22383029285599393, 2.01, 0.3);

            
        testMultiple(1, 0.2, 2, 0.1, 
            2, 0.4127953488110059, 2.02, 0.5);
        testMultiple(1, 0.2, -2, 0.1, 
            -2, 0.4127953488110059, -2.1, 0.42);
        testMultiple(-1, 0.2, 2, 0.1, 
            -2, 0.4127953488110059, -1.9, 0.38);
        testMultiple(-1, 0.2, -2, 0.1, 
            2, 0.4127953488110059, 2.02, 0.5);

            
        testMultiple(1, 0.2, 0, 0.1, 
            0, 0.10198039027185571, 0, 0.12);
        testMultiple(-1, 0.2, 0, 0.1, 
            0, 0.10198039027185571, 0.01999999999999999, 0.1);

        
    }


    @Test
    public void testLargeDiff() {
        assertEquals(13316075197586562L, 64919121L*205117922L);
        assertEquals(13316075197586561L, 159018721L*83739041L);
        assertEquals(1L, 64919121L*205117922L - 159018721L*83739041L);

        testMultiple(64919121, Double.NaN, 205117922, Double.NaN, 
                    6658037598793281L*2, 2.8284271247461903, 6658037598793281L*2, 4); 
        testMultiple(-159018721L, Double.NaN, 83739041, Double.NaN, 
                    -6658037598793280L*2, 3.4641016151377544, -6658037598793280L*2, 3);

        try {
            final VarDbl op1 = new VarDbl( 6658037598793281L*2, 2.8284271247461903*2.8284271247461903);
            final VarDbl op2 = new VarDbl(-6658037598793280L*2, 3.4641016151377544*3.4641016151377544);
            final VarDbl res = (VarDbl) op1.add(op2);
            assertEquals(2, res.value(), TOLERANCE);
            assertEquals(4.472135954999579393, res.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException | TypeException e) {
            fail();
        }
        
        try {
            final IntvDbl op1 = new IntvDbl( 6658037598793281L*2, 4);
            final IntvDbl op2 = new IntvDbl(-6658037598793280L*2, 3);
            final IntvDbl res = (IntvDbl) op1.add(op2);
            assertEquals(2, res.value(), TOLERANCE);
            assertEquals(7, res.uncertainty(), TOLERANCE);
        } catch (ValueException | UncertaintyException | TypeException e) {
            fail();
        }

    }
    
}
