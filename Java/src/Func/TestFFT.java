package Func;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TestFFT {
    static final double TOLERANCE = 2E-16;
    
    static final double q1 = Math.sin(Math.PI/8);
    static final double q2 = Math.sin(Math.PI/4);
    static final double q3 = Math.sin(Math.PI/8*3);

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
}
