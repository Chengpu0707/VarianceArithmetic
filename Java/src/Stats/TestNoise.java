package Stats;

import org.apache.commons.math3.distribution.NormalDistribution;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class TestNoise {
    final double TOLERANCE = 3E-16;
    final Noise noise = new Noise();
    private static final NormalDistribution distr = new NormalDistribution();

    @Test
    public void testGaussan_dev1() {
        // demonstrate the stability of gaussian noise using 100,000 samples
        final Histogram histo = new Histogram(3, 5);
        final double[] sExpected = new double[31];
        for (int i = 0; i < 31; ++i) {
            final double x0 = -3.0 + i * 0.2;
            sExpected[i] = distr.cumulativeProbability(x0 + 0.2) - distr.cumulativeProbability(x0);
        }
        for (int i = 0; i < 100000; ++i) {
            histo.accum(noise.gaussian(1), i + 1);
        }
        assertEquals(0, histo.stat().avg(), 3E-2); 
        assertEquals(1, histo.stat().dev(), 5E-2); 
        assertTrue(0 < histo.lower());
        assertTrue(0 < histo.upper());
        final double[] sHisto = histo.histo();
        for (int i = 10; i < sExpected.length - 10; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.2);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
        for (int i = 5; i < 10; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.3);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
        for (int i = sExpected.length - 10; i < sExpected.length - 5; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.3);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
    }

    @Test
    public void testGaussan_dev2() {
        // demonstrate the stability of gaussian noise using 100,000 samples
        final Histogram histo = new Histogram(3, 5);
        final double[] sExpected = new double[31];
        for (int i = 0; i < 31; ++i) {
            final double x0 = -1.5 + i * 0.1;
            sExpected[i] = distr.cumulativeProbability(x0 + 0.1) - distr.cumulativeProbability(x0);
        }
        for (int i = 0; i < 200000; ++i) {
            histo.accum(noise.gaussian(2), i + 1);
        }
        assertEquals(0, histo.stat().avg(), 3E-2); 
        assertEquals(2, histo.stat().dev(), 5E-2); 
        assertTrue(0 < histo.lower());
        assertTrue(0 < histo.upper());
        final double[] sHisto = histo.histo();
        for (int i = 10; i < sExpected.length - 10; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.2);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
        for (int i = 5; i < 10; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.3);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
        for (int i = sExpected.length - 10; i < sExpected.length - 5; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.3);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
    }

    @Test
    public void testWite_dev1() {
        final Histogram histo = new Histogram(2, 5);
        for (int i = 0; i < 100000; ++i) {
            histo.accum(noise.white(1), i);
        }
        final double[] sExpected = new double[21];
        for (int i = 2; i < 19; ++i) {
            sExpected[i] = 1 / Math.sqrt(3) / 10;
        }
        assertTrue(100000 == histo.stat().count());
        assertEquals(0, histo.stat().avg(), 5E-2); 
        assertEquals(1, histo.stat().dev(), 3E-2); 
        assertEquals(0, histo.lower(), TOLERANCE);
        assertEquals(0, histo.upper(), TOLERANCE);
        assertEquals(Math.sqrt(3), -histo.stat().min(), 5E-2);
        assertEquals(Math.sqrt(3), histo.stat().max(), 5E-2);
        final double[] sHisto = histo.histo();
        for (int i = 2; i < sExpected.length - 2; ++i) {  
            try {
                assertEquals(1, sHisto[i] / sExpected[i], 0.05);
            } catch (AssertionError e) {
                System.out.printf("At index %d, %f vs %f", i, sHisto[i], sExpected[i]);
                throw e;
            }  
        }
    }
}
