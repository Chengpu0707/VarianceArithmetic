package Func;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TestFFT {
    static final double TOLERANCE = 2E-16;
    
    @Test
    public void testSine() {
        final double q1 = Math.sin(Math.PI/8);
        final double q2 = Math.sin(Math.PI/4);
        final double q3 = Math.sin(Math.PI/8*3);

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
}
