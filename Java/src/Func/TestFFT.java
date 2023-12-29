package Func;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal;
import Type.VarDbl;
import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class TestFFT {
    static final double TOLERANCE = 2E-16;
    
    static final double q1 = Math.sin(Math.PI/8);
    static final double q2 = Math.sin(Math.PI/4);
    static final double q3 = Math.sin(Math.PI/8*3);

    static final VarDbl ZERO = new VarDbl();
    static final VarDbl POS_ONE = new VarDbl(1);
    static final VarDbl NEG_ONE = new VarDbl(-1);

    @Test
    public void testSine() {
        assertEquals(-q3, FFT.sin(-5, 4), TOLERANCE);

        assertEquals(-1,  FFT.sin(-4, 4), TOLERANCE);
        assertEquals(-q3, FFT.sin(-3, 4), TOLERANCE);
        assertEquals(-q2, FFT.sin(-2, 4), TOLERANCE);
        assertEquals(-q1, FFT.sin(-1, 4), TOLERANCE);

        assertEquals(0, FFT.sin(0, 4), TOLERANCE);
        assertEquals(q1, FFT.sin(1, 4), TOLERANCE);
        assertEquals(q2, FFT.sin(2, 4), TOLERANCE);
        assertEquals(q3, FFT.sin(3, 4), TOLERANCE);

        assertEquals(1, FFT.sin(4, 4), TOLERANCE);
        assertEquals(q3, FFT.sin(5, 4), TOLERANCE);
        assertEquals(q2, FFT.sin(6, 4), TOLERANCE);
        assertEquals(q1, FFT.sin(7, 4), TOLERANCE);

        assertEquals(0, FFT.sin(8, 4), TOLERANCE);
        assertEquals(-q1, FFT.sin(9, 4), TOLERANCE);
        assertEquals(-q2, FFT.sin(10, 4), TOLERANCE);
        assertEquals(-q3, FFT.sin(11, 4), TOLERANCE);
        assertEquals(-1, FFT.sin(12, 4), TOLERANCE);

        assertEquals(-q3, FFT.sin(13, 4), TOLERANCE);
        assertEquals(-q2, FFT.sin(14, 4), TOLERANCE);
        assertEquals(-q1, FFT.sin(15, 4), TOLERANCE);
        assertEquals(0, FFT.sin(16, 4), TOLERANCE);

        assertEquals(q1, FFT.sin(17, 4), TOLERANCE);
        assertEquals(q2, FFT.sin(18, 4), TOLERANCE);
        assertEquals(q3, FFT.sin(19, 4), TOLERANCE);
        
        assertEquals(1, FFT.sin(20, 4), TOLERANCE);
    }

    @Test
    public void testCosine() {
        assertEquals(-q1, FFT.cos(-5, 4), TOLERANCE);

        assertEquals(0, FFT.cos(-4, 4), TOLERANCE);
        assertEquals(q1, FFT.cos(-3, 4), TOLERANCE);
        assertEquals(q2, FFT.cos(-2, 4), TOLERANCE);
        assertEquals(q3, FFT.cos(-1, 4), TOLERANCE);

        assertEquals(1, FFT.cos(0, 4), TOLERANCE);
        assertEquals(q3, FFT.cos(1, 4), TOLERANCE);
        assertEquals(q2, FFT.cos(2, 4), TOLERANCE);
        assertEquals(q1, FFT.cos(3, 4), TOLERANCE);

        assertEquals(0, FFT.cos(4, 4), TOLERANCE);
        assertEquals(-q1, FFT.cos(5, 4), TOLERANCE);
        assertEquals(-q2, FFT.cos(6, 4), TOLERANCE);
        assertEquals(-q3, FFT.cos(7, 4), TOLERANCE);

        assertEquals(-1, FFT.cos(8, 4), TOLERANCE);
        assertEquals(-q3, FFT.cos(9, 4), TOLERANCE);
        assertEquals(-q2, FFT.cos(10, 4), TOLERANCE);
        assertEquals(-q1, FFT.cos(11, 4), TOLERANCE);

        assertEquals(0, FFT.cos(12, 4), TOLERANCE);
        assertEquals(q1, FFT.cos(13, 4), TOLERANCE);
        assertEquals(q2, FFT.cos(14, 4), TOLERANCE);
        assertEquals(q3, FFT.cos(15, 4), TOLERANCE);

        assertEquals(1, FFT.cos(16, 4), TOLERANCE);
        assertEquals(q3, FFT.cos(17, 4), TOLERANCE);
        assertEquals(q2, FFT.cos(18, 4), TOLERANCE);
        assertEquals(q1, FFT.cos(19, 4), TOLERANCE);

        assertEquals(0, FFT.cos(20, 4), TOLERANCE);
    }

    @Test
    public void testBitReverseIndices() {
        assertArrayEquals( new int[] {0, 2, 1, 3}, FFT.bitReversedIndices(2));
        assertArrayEquals( new int[] {0, 4, 2, 6, 1, 5, 3, 7}, FFT.bitReversedIndices(3));
        assertArrayEquals( new int[] {0, 8, 4, 12, 2, 10, 6, 14, 
                                      1, 9, 5, 13, 3, 11, 7, 15}, FFT.bitReversedIndices(4));
                                      
        for (int order = 2; order <= FFT.MAX_ORDER; ++order) {
            final int[] sRes = FFT.bitReversedIndices(order);
            assertNotNull(String.format("order=%d", order), sRes);
            for (int i = 0; i < sRes.length; ++i) {
                final int j = sRes[i];
                int br = 0, org = i;
                for (int k = 0; k < order; ++k) {
                    br <<= 1;
                    br |= (org & 1);
                    org >>= 1;
                }
                assertEquals(String.format("order=%d, index=%d: %d != %d", order, i, br, j), br, j);
            }
        }
    }

    void testSin(int order, int freq, double sigma) {
        final int sz = 1 << order;
        assertTrue(String.format("Sin invalid freq %d for order %d", freq, order), 
            freq * 2 < sz);
        final VarDbl[] sData = new VarDbl[sz << 1];
        for (int i = 0; i < sz; ++i) {
            sData[i << 1] = new VarDbl( FFT.sin(freq * i, order) );
            sData[(i << 1) + 1] = ZERO;
        }
        try {
            final IReal[] sSpec = FFT.transform(sData, true);
            for (int i = 0; i < sz; ++i) {
                if (i == ((freq << 1) + 1))
                    assertEquals(String.format("Sin Forward FFT order %d freq %d, index %d: %d != %s", 
                            order, freq, i, sz >> 1, sSpec[i]),
                        sz >> 1, sSpec[i].value(), 
                        Math.abs(sigma * sSpec[i].value() * sSpec[i].uncertainty()));
                else if (i == (((sz - freq) << 1) + 1))
                    assertEquals(String.format("Sin Forward FFT order %d freq %d, index %d: %d != %s", 
                            order, freq, i, -(sz >> 1), sSpec[i]),
                        -(sz >> 1), sSpec[i].value(), 
                        Math.abs(sigma * sSpec[i].value() * sSpec[i].uncertainty()));
                else
                    assertEquals(String.format("Sin Forward FFT order %d freq %d, index %d: %d != %s", 
                            order, freq, i, sz >> 1, sSpec[i]),
                        0, sSpec[i].value(), sigma * sSpec[i].uncertainty());
             }

            final IReal[] sRev = FFT.transform(sSpec, false);
            for (int i = 0; i < 8; ++i) {
                assertEquals(String.format("Sin Reverse FFT order %d freq %d, index %d: %s != %s", 
                        order, freq, i, sRev[i].toString(), sData[i].toString()),
                    sData[i].value(), sRev[i].value(), 
                    (sRev[i].value() == 0)
                        ? sigma * sRev[i].uncertainty()
                        : Math.abs(sigma * sRev[i].value() * sRev[i].uncertainty()));
            }
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    void testCos(int order, int freq, double sigma) {
        final int sz = 1 << order;
        assertTrue(String.format("Cos invalid freq %d for order %d", freq, order), 
            freq * 2 < sz);
        final VarDbl[] sData = new VarDbl[sz << 1];
        for (int i = 0; i < sz; ++i) {
            sData[i << 1] = new VarDbl( FFT.cos(freq * i, order) );
            sData[(i << 1) + 1] = ZERO;
        }
        try {
            final IReal[] sSpec = FFT.transform(sData, true);
            for (int i = 0; i < sz; ++i) {
                if ((i == (freq << 1)) || (i == ((sz - freq) << 1)))
                    assertEquals(String.format("Cos Forward FFT order %d freq %d, index %d: %d != %s", 
                            order, freq, i, sz >> 1, sSpec[i]),
                        sz >> 1, sSpec[i].value(), 
                        Math.abs(sigma * sSpec[i].value() * sSpec[i].uncertainty()));
                else
                    assertEquals(String.format("Sin Forward FFT order %d freq %d, index %d: %d != %s", 
                            order, freq, i, sz >> 1, sSpec[i]),
                        0, sSpec[i].value(), sSpec[i].uncertainty());
             }

            final IReal[] sRev = FFT.transform(sSpec, false);
            for (int i = 0; i < 8; ++i) {
                assertEquals(String.format("Cos Reverse FFT order %d freq %d, index %d: %s != %s", 
                        order, freq, i, sRev[i].toString(), sData[i].toString()),
                    sData[i].value(), sRev[i].value(), 
                    (sRev[i].value() == 0)
                        ? sigma * sRev[i].uncertainty()
                        : Math.abs(sigma * sRev[i].value() * sRev[i].uncertainty()));
            }
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testAllTransform() {
        final double sigma = 3.0;
        testCos(3, 1, sigma);
        for (int order = 2; order <= FFT.MAX_ORDER; ++order) {
            final int maxFreq = 1 << (order - 1);
            for (int freq = 1; freq < maxFreq; ++ freq) {
                testSin(order, freq, sigma);
                testCos(order, freq, sigma);
            }
        }
    }
}
