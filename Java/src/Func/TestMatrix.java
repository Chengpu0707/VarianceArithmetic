package Func;

/**
 * JUnit tests for Matrix — construction, get/set, bottom-up Laplace
 * determinant, adjugate (Formula 6.1) and its variance (Formula 6.9),
 * the M·adj(M) = |M|·I identity, variance propagation, the random-matrix
 * factory, and a Monte-Carlo diagnostic confirming Formula 6.9 predicts
 * the unconditional adjugate-cell variance correctly.
 */
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.junit.Test;

import Stats.Stat;
import Type.InitException;
import Type.VarDbl;


public class TestMatrix {

    private static void put(Matrix m, int r, int c, double val) throws InitException {
        m.set(r, c, new VarDbl(val, 0));
    }

    @Test
    public void testConstruction() throws InitException {
        final Matrix m = new Matrix(3);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                assertEquals(0.0, m.get(r, c).value(),       0);
                assertEquals(0.0, m.get(r, c).uncertainty(), 0);
            }
    }

    @Test
    public void testGetSet() throws InitException {
        final Matrix m = new Matrix(2);
        put(m, 0, 0, 1.0); put(m, 0, 1, 2.0);
        put(m, 1, 0, 3.0); put(m, 1, 1, 4.0);
        assertEquals(1.0, m.get(0, 0).value(), 0);
        assertEquals(2.0, m.get(0, 1).value(), 0);
        assertEquals(3.0, m.get(1, 0).value(), 0);
        assertEquals(4.0, m.get(1, 1).value(), 0);
    }

    @Test
    public void testDetermOne() throws InitException {
        final Matrix m = new Matrix(1);
        put(m, 0, 0, 5.0);
        assertEquals(5.0, m.determ().value(), 1e-12);
    }

    @Test
    public void testDetermIdentity() throws InitException {
        final Matrix m = new Matrix(4);
        for (int i = 0; i < 4; ++i) put(m, i, i, 1.0);
        assertEquals(1.0, m.determ().value(), 1e-12);
    }

    @Test
    public void testDetermZero() throws InitException {
        assertEquals(0.0, new Matrix(3).determ().value(), 1e-12);
    }

    @Test
    public void testDetermDiagonal() throws InitException {
        final Matrix m = new Matrix(4);
        put(m, 0, 0, 2.0); put(m, 1, 1, 3.0);
        put(m, 2, 2, 5.0); put(m, 3, 3, 7.0);
        assertEquals(2.0 * 3.0 * 5.0 * 7.0, m.determ().value(), 1e-9);
    }

    @Test
    public void testDeterm2x2() throws InitException {
        final Matrix m = new Matrix(2);
        put(m, 0, 0, 1.0); put(m, 0, 1, 2.0);
        put(m, 1, 0, 3.0); put(m, 1, 1, 4.0);
        assertEquals(-2.0, m.determ().value(), 1e-12);   // 1·4 − 2·3
    }

    @Test
    public void testDeterm3x3() throws InitException {
        // Expand along row 0: 6·(-38) − 1·(-14) + 36 = −178.
        final Matrix m = new Matrix(3);
        put(m, 0, 0,  6.0); put(m, 0, 1,  1.0); put(m, 0, 2,  1.0);
        put(m, 1, 0,  4.0); put(m, 1, 1, -2.0); put(m, 1, 2,  5.0);
        put(m, 2, 0,  2.0); put(m, 2, 1,  8.0); put(m, 2, 2, -1.0);
        assertEquals(-178.0, m.determ().value(), 1e-9);
    }

    @Test
    public void testDeterm4x4UpperTri() throws InitException {
        final Matrix m = new Matrix(4);
        put(m, 0, 0, 1.0); put(m, 0, 1, 5.0); put(m, 0, 2, 9.0); put(m, 0, 3, 2.0);
                           put(m, 1, 1, 2.0); put(m, 1, 2, 6.0); put(m, 1, 3, 7.0);
                                              put(m, 2, 2, 3.0); put(m, 2, 3, 8.0);
                                                                 put(m, 3, 3, 4.0);
        assertEquals(24.0, m.determ().value(), 1e-9);
    }

    @Test
    public void testAdjugate2x2() throws InitException {
        // adj([[a,b],[c,d]]) = [[d,-b],[-c,a]]
        final Matrix m = new Matrix(2);
        put(m, 0, 0, 1.0); put(m, 0, 1, 2.0);
        put(m, 1, 0, 3.0); put(m, 1, 1, 4.0);
        final Matrix adj = m.adjugate();
        assertEquals( 4.0, adj.get(0, 0).value(), 1e-12);
        assertEquals(-2.0, adj.get(0, 1).value(), 1e-12);
        assertEquals(-3.0, adj.get(1, 0).value(), 1e-12);
        assertEquals( 1.0, adj.get(1, 1).value(), 1e-12);
    }

    @Test
    public void testAdjugateIdentity() throws InitException {
        final Matrix m = new Matrix(3);
        for (int i = 0; i < 3; ++i) put(m, i, i, 1.0);
        final Matrix adj = m.adjugate();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                assertEquals(i == j ? 1.0 : 0.0, adj.get(i, j).value(), 1e-12);
    }

    @Test
    public void testAdjugateRoundtrip3x3() throws InitException {
        // Verify M · adj(M) = |M| · I on a non-trivial 3×3.
        final Matrix m = new Matrix(3);
        put(m, 0, 0,  6.0); put(m, 0, 1,  1.0); put(m, 0, 2,  1.0);
        put(m, 1, 0,  4.0); put(m, 1, 1, -2.0); put(m, 1, 2,  5.0);
        put(m, 2, 0,  2.0); put(m, 2, 1,  8.0); put(m, 2, 2, -1.0);
        final double det = m.determ().value();
        final Matrix adj = m.adjugate();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double acc = 0;
                for (int k = 0; k < 3; ++k) acc += m.get(i, k).value() * adj.get(k, j).value();
                assertEquals(i == j ? det : 0.0, acc, 1e-10);
            }
    }

    @Test
    public void testDetermPropagatesVariance() throws InitException {
        // 2×2 diagonal of independent uncertain entries — determinant gets non-zero unc.
        final Matrix m = new Matrix(2);
        m.set(0, 0, new VarDbl(3.0, 0.1));
        m.set(1, 1, new VarDbl(5.0, 0.1));
        final VarDbl d = m.determ();
        assertEquals(15.0, d.value(), 1e-12);
        assertTrue("uncertainty should propagate", d.uncertainty() > 0);
    }

    @Test
    public void testRandomMatrix() throws InitException {
        final double dev = 1e-3;
        final Matrix m = Matrix.randomMatrix(4, dev);
        final double expectedUnc = 256.0 * dev;
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                assertEquals(expectedUnc, m.get(r, c).uncertainty(), 1e-12);
                final double v = m.get(r, c).value();
                assertTrue("cell value out of expected range: " + v, v > -1024 && v < 1024);
            }
        }
        assertNotNull(m.determ());
    }

    /**
     * Confirms Formula 6.9 implementation is faithful: the unconditional
     * predicted unc (from a precise-int matrix with uncertainty `dev`)
     * matches the empirical std of (adj_noisy − adj_precise) across samples.
     * Normalized error std should sit at ≈ 1 at every noise level.
     */
    @Test
    public void testAdjugateUncertaintyVsMonteCarlo() throws InitException {
        final int[][] sInt = {{2, -1, 3}, {0, 4, -2}, {1, -3, 1}};
        final int size = sInt.length;
        final int samples = 1000;
        final Matrix adjPrecise = Matrix.asMatrix(sInt).adjugate();
        final double[] noises = {1e-3, 1e-2, 1e-1, 1.0, 10.0};

        System.out.println();
        System.out.printf("Monte-Carlo adjugate diagnostic (size %d):%n", size);
        System.out.printf("%8s  %12s  %14s  %14s  %14s%n",
                          "noise", "dev", "MC_norm_std", "MC_val_std", "pred_unc");

        for (final double noiseLevel : noises) {
            final double dev = 256.0 / Math.sqrt(3.0) * noiseLevel;

            // Predictor: precise int means with uncertainty `dev` — Formula 6.9
            // applied here gives the true unconditional variance of adj(M).
            final Matrix preciseWithDev = new Matrix(size);
            for (int r = 0; r < size; ++r)
                for (int c = 0; c < size; ++c)
                    preciseWithDev.set(r, c, new VarDbl((double) sInt[r][c], dev));
            final Matrix adjPredicted = preciseWithDev.adjugate();

            final Random rng = new Random(42);
            final Stat normErr = new Stat();
            final Stat valErr  = new Stat();
            for (int s = 0; s < samples; ++s) {
                final Matrix adjNoisy = Matrix.addNoise(sInt, dev, rng).adjugate();
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < size; ++j) {
                        final double err = adjNoisy.get(i, j).value()
                                          - adjPrecise.get(i, j).value();
                        final double unc = adjPredicted.get(i, j).uncertainty();
                        valErr.accum(err, 0);
                        if (unc > 0) normErr.accum(err / unc, 0);
                    }
                }
            }

            double sumUnc = 0;
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j)
                    sumUnc += adjPredicted.get(i, j).uncertainty();
            System.out.printf("%8.0e  %12.3e  %14.4f  %14.4e  %14.4e%n",
                              noiseLevel, dev, normErr.dev(), valErr.dev(), sumUnc / (size * size));

            final double r = normErr.dev();
            assertTrue("normalized std at noise=" + noiseLevel + " should be ≈ 1, was " + r,
                       r > 0.85 && r < 1.15);
        }
    }
}
