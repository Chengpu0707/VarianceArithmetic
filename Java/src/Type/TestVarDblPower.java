package Type;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

/*
public class TestVarDblPower extends VarDbl {
    final static double TOLERANCE = 3E-16;

    @Test
    public void testTaylor() {
        double[] sTaylor;

        sTaylor = powerTaylor(1);
        assertEquals(1, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);

        sTaylor = powerTaylor(2);
        assertEquals(2, sTaylor[1], TOLERANCE);
        assertEquals(1, sTaylor[2], TOLERANCE);
        assertEquals(0, sTaylor[3], TOLERANCE);

        sTaylor = powerTaylor(0.5);
        assertEquals(1.0/2, sTaylor[1], TOLERANCE);
        assertEquals(-1.0/8, sTaylor[2], TOLERANCE);
        assertEquals(1.0/16, sTaylor[3], TOLERANCE);
        assertEquals(-5.0/128, sTaylor[4], TOLERANCE);

        sTaylor = powerTaylor(-1);
        for (int i = 1; i < sTaylor.length; ++i) {
            assertEquals(((i % 2) == 1)? -1 : 1, sTaylor[i], TOLERANCE);
        }
    }

    @Test
    public void testGetPowerCoeff() {
        double[] sCoeff; 
        
        sCoeff = getPowerCoeff(1);
        assertEquals(0, sCoeff[1], TOLERANCE);
        assertEquals(1, sCoeff[2], TOLERANCE);
        assertEquals(0, sCoeff[3], TOLERANCE);

        sCoeff = getPowerCoeff(2);
        assertEquals(0, sCoeff[1], TOLERANCE);
        assertEquals(4, sCoeff[2], TOLERANCE);
        assertEquals(4, sCoeff[3], TOLERANCE);
        assertEquals(1, sCoeff[4], TOLERANCE);
        assertEquals(0, sCoeff[5], TOLERANCE);

        sCoeff = getPowerCoeff(0.5);
        assertEquals(0, sCoeff[1], TOLERANCE);
        assertEquals(1.0/4, sCoeff[2], TOLERANCE);
        assertEquals(-1.0/8, sCoeff[3], TOLERANCE);
        assertEquals(1.0/16 + 1.0/64, sCoeff[4], TOLERANCE);
        assertEquals(-5.0/128 - 1.0/64, sCoeff[5], TOLERANCE);

        sCoeff = getPowerCoeff(-1);
        for (int i = 1; i < sCoeff.length; ++i) {
            assertEquals(((i % 2) == 1)? -(i - 1) : (i - 1), sCoeff[i], TOLERANCE);
        }
    }
    
}

*/
