package Stats;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class TestNoise {
    final double TOLERANCE = 3E-16;
    final Noise noise = new Noise();

    @Test
    public void testGaussan_dev0() {
        final Histogram histo = new Histogram(3, 4);
        for (int i = 0; i < 1000; ++i) {
            histo.accum(noise.gaussian(0));
        }
        assertEquals(1000, histo.stat().count());
        assertEquals(0, histo.stat().avg(), TOLERANCE); 
        assertEquals(0, histo.stat().dev(), TOLERANCE); 
        assertEquals(0, histo.stat().min(), TOLERANCE);
        assertEquals(0, histo.stat().max(), TOLERANCE);
        assertEquals(0, histo.actRange(), TOLERANCE);
    }

    @Test
    public void testGaussan_dev1() {
        final Histogram histo = new Histogram(3, 4);
        for (int i = 0; i < 10000; ++i) {
            histo.accum(noise.gaussian(1));
        }
        assertTrue(10000 >= histo.stat().count());
        assertTrue(9900 <= histo.stat().count());
        assertEquals(0, histo.stat().avg(), 3E-2); 
        assertEquals(1, histo.stat().dev(), 5E-2); 
        final int actRange = (int) Math.round(Math.max( -histo.stat().min(), histo.stat().max()) * histo.divids());
        assertEquals(actRange, histo.actRange());
        assertArrayEquals(new double[]{0.06, 0.075, 0.089, 0.098, 0.10, 0.098, 0.089, 0.075, 0.06}, 
                histo.histo(4), 1E-2);
    }

    @Test
    public void testGaussan_dev2() {
        final Histogram histo = new Histogram(6, 4);
        for (int i = 0; i < 20000; ++i) {
            histo.accum(noise.gaussian(2));
        }
        assertTrue(20000 >= histo.stat().count());
        assertTrue(19900 <= histo.stat().count());
        assertEquals(0, histo.stat().avg(), 5E-2); 
        assertEquals(2, histo.stat().dev(), 5E-2); 
        final int actRange = (int) Math.round(Math.max( -histo.stat().min(), histo.stat().max()) * histo.divids());
        assertEquals(actRange, histo.actRange());
        assertArrayEquals(new double[]{0.043, 0.046, 0.049, 0.049, 0.052, 0.049, 0.049, 0.046, 0.043}, 
                histo.histo(4), 1E-2);
    }

    @Test
    public void testWite_dev2() {
        final Histogram histo = new Histogram(4, 2);
        for (int i = 0; i < 20000; ++i) {
            histo.accum(noise.white(2));
        }
        assertTrue(20000 == histo.stat().count());
        assertEquals(0, histo.stat().avg(), 5E-2); 
        assertEquals(2, histo.stat().dev(), 3E-2); 
        final int actRange = (int) Math.round(Math.max( -histo.stat().min(), histo.stat().max()) * histo.divids());
        assertEquals(actRange, histo.actRange());
        assertEquals(Math.sqrt(3)*2, -histo.stat().min(), 5E-2);
        assertEquals(Math.sqrt(3)*2, histo.stat().max(), 5E-2);
        assertArrayEquals(new double[]{0.071, 0.071, 0.071, 0.071, 0.071, 0.071, 0.071, 0.071, 0.071}, 
                histo.histo(4), 1E-2);
    }
}
