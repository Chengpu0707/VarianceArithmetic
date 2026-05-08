package Func;

/**
 * JUnit tests for {@link Matrix} — verifies permutation sign, square-matrix
 * checks, multiply, linear, and adjugate / adjugate_mul on Long, Double,
 * Fraction, and VarDbl element types. Mirrors a deterministic subset of
 * Python testMatrix.py.
 *
 * Intentionally omitted (require Python-specific RNG/numpy plumbing or are
 * Monte-Carlo coverage runs):
 *   - testPermut: checks Python itertools order, not relevant to Java.
 *   - testNumpyIntMatrix, testHilbert: depend on numpy.
 *   - testCreateIntMatrix, testAddNoise: live RNG statistical thresholds.
 *   - testFractions, testInt (TestFractions): trivial Python operator tests.
 *   - Adjugate.testAdjugate, testIdealCoverage, testProperCoverage: slow
 *     Monte-Carlo dump tests.
 *   - testException, testInversionException: rely on division/overflow paths
 *     in VarDbl that mirror Python-specific exceptions.
 */
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.apache.commons.math3.fraction.Fraction;
import org.junit.Test;

import Type.InitException;
import Type.VarDbl;


public class TestMatrix {

    private static final double FLOAT_EPS = 1e-7;

    /* ----------------------- helpers ----------------------- */

    private static Object[][] mat(Object... row) {
        // Build a square matrix from a flat list. row.length must be a perfect square.
        int total = row.length;
        int n = (int) Math.round(Math.sqrt(total));
        if (n * n != total)
            throw new IllegalArgumentException("Not a square count: " + total);
        Object[][] out = new Object[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                out[i][j] = row[i * n + j];
        return out;
    }

    private static Long L(long v) { return Long.valueOf(v); }

    private void verifyValue(Object expected, Object actual) {
        if (expected instanceof VarDbl) {
            assertTrue("Expected VarDbl, got " + actual, actual instanceof VarDbl);
            VarDbl ve = (VarDbl) expected;
            VarDbl va = (VarDbl) actual;
            assertEquals(ve.value(), va.value(), Math.max(1e-7, Math.abs(ve.value()) * 1e-7));
            double evar = ve.variance();
            double avar = va.variance();
            assertEquals(evar, avar, Math.max(1e-12, Math.abs(evar) * 1e-2));
        } else if (expected instanceof Double || expected instanceof Float) {
            double e = ((Number) expected).doubleValue();
            double a;
            if (actual instanceof Number) a = ((Number) actual).doubleValue();
            else if (actual instanceof Fraction) a = ((Fraction) actual).doubleValue();
            else { fail("Cannot compare " + expected + " to " + actual); return; }
            assertEquals(e, a, Math.max(FLOAT_EPS, Math.abs(e) * 1e-7));
        } else if (expected instanceof Fraction && actual instanceof Fraction) {
            assertEquals(expected, actual);
        } else if (expected instanceof Long && actual instanceof Long) {
            assertEquals(((Long) expected).longValue(), ((Long) actual).longValue());
        } else if (expected instanceof Long && actual instanceof Fraction) {
            assertEquals(new Fraction(((Long) expected).intValue(), 1), actual);
        } else {
            assertEquals(expected, actual);
        }
    }

    private void verifyIdentity(Object det, Object[][] ssId) {
        verifyIdentity(det, ssId, 1e-6, 3, 5);
    }

    private void verifyIdentity(Object det, Object[][] ssId,
                                double delta, int uncRange, int places) {
        int size = ssId.length;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Object cell = ssId[i][j];
                if (cell instanceof VarDbl) {
                    double dv = (det instanceof VarDbl) ? ((VarDbl) det).value()
                            : ((Number) det).doubleValue();
                    double expected = (i == j) ? dv : 0.0;
                    double actual = ((VarDbl) cell).value();
                    double tol = Math.max(Math.abs(dv) * delta, ((VarDbl) cell).uncertainty() * uncRange);
                    if (tol == 0) tol = 1e-9;
                    assertEquals(String.format("ssId[%d][%d]", i, j), expected, actual, tol);
                } else if (cell instanceof Double) {
                    double dv = ((Number) det).doubleValue();
                    double expected = (i == j) ? dv : 0.0;
                    double actual = ((Double) cell).doubleValue();
                    assertEquals(String.format("ssId[%d][%d]", i, j),
                            expected, actual, Math.pow(10, -places));
                } else if (cell instanceof Fraction) {
                    Fraction expected = (i == j)
                        ? ((det instanceof Fraction) ? (Fraction) det : Matrix.objectEquals(det, L(0)) ? Fraction.ZERO : new Fraction(((Number) det).intValue(), 1))
                        : Fraction.ZERO;
                    assertEquals(expected, cell);
                } else {
                    long expected = (i == j) ? ((Number) det).longValue() : 0L;
                    long actual = ((Number) cell).longValue();
                    assertEquals(String.format("ssId[%d][%d]", i, j), expected, actual);
                }
            }
        }
    }

    private void verify(Object[][] ssMat, Object det, Object[][] ssAdj) throws InitException {
        Object[][] prod = Matrix.multiply(ssMat, ssAdj);
        verifyIdentity(det, prod);
        Matrix.AdjugateResult r = Matrix.adjugate(ssMat);
        verifyValue(det, r.det);
        verifyIdentity(det, Matrix.multiply(ssMat, r.adj));
        int size = ssMat.length;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                verifyValue(ssAdj[i][j], r.adj[i][j]);
            }
        }
    }


    /* ----------------------- TestPermutation ----------------------- */

    @Test
    public void testSign() {
        assertEquals(1, Matrix.permutSign(new int[] {0, 1}));
        assertEquals(-1, Matrix.permutSign(new int[] {1, 0}));

        assertEquals(1, Matrix.permutSign(new int[] {0, 1, 2}));
        assertEquals(-1, Matrix.permutSign(new int[] {0, 2, 1}));
        assertEquals(-1, Matrix.permutSign(new int[] {1, 0, 2}));
        assertEquals(1, Matrix.permutSign(new int[] {1, 2, 0}));
        assertEquals(1, Matrix.permutSign(new int[] {2, 0, 1}));
        assertEquals(-1, Matrix.permutSign(new int[] {2, 1, 0}));
    }

    @Test
    public void testPermutationCount() {
        assertEquals(2, Matrix.permutations(2).size());
        assertEquals(6, Matrix.permutations(3).size());
        assertEquals(24, Matrix.permutations(4).size());
    }

    @Test
    public void testCombinationCount() {
        assertEquals(3, Matrix.combinations(3, 1).size());
        assertEquals(3, Matrix.combinations(3, 2).size());
        assertEquals(1, Matrix.combinations(3, 3).size());
    }


    /* ----------------------- TestSquareMatrix ----------------------- */

    @Test
    public void testIsSquareMatrix() throws InitException {
        Object[][] ok = new Object[][] {
            { new VarDbl(1.0), Double.valueOf(2.0) },
            { L(3), new VarDbl(4.0, 1.0) },
        };
        assertTrue(Matrix.isSquareMatrix(ok));

        // null element rejected.
        Object[][] nullEntry = new Object[][] {
            { new VarDbl(1.0), Double.valueOf(2.0) },
            { L(3), null },
        };
        assertFalse(Matrix.isSquareMatrix(nullEntry));

        // restricted type set rejects VarDbl.
        Class<?>[] intFloat = new Class<?>[] { Long.class, Integer.class, Double.class, Float.class };
        assertFalse(Matrix.isSquareMatrix(ok, intFloat));

        // non-square rejected.
        Object[][] nonSquare = new Object[][] {
            { L(1), L(2) },
            { L(3) },
        };
        assertFalse(Matrix.isSquareMatrix(nonSquare));
    }


    /* ----------------------- TestMultiply ----------------------- */

    @Test
    public void testIntSize2() throws InitException {
        Object[][] ssId = mat(L(1), L(0), L(0), L(1));
        Object[][] ssTr = mat(L(0), L(1), L(1), L(0));
        Object[][] ssMat = mat(L(1), L(-2), L(-3), L(4));

        assertTrue(Matrix.elementsEqual(ssId, Matrix.multiply(ssId, ssId)));

        assertTrue(Matrix.elementsEqual(ssTr, Matrix.multiply(ssTr, ssId)));
        assertTrue(Matrix.elementsEqual(ssTr, Matrix.multiply(ssId, ssTr)));
        assertTrue(Matrix.elementsEqual(ssId, Matrix.multiply(ssTr, ssTr)));

        assertTrue(Matrix.elementsEqual(ssMat, Matrix.multiply(ssId, ssMat)));
        assertTrue(Matrix.elementsEqual(ssMat, Matrix.multiply(ssMat, ssId)));

        assertTrue(Matrix.elementsEqual(mat(L(-3), L(4), L(1), L(-2)), Matrix.multiply(ssTr, ssMat)));
        assertTrue(Matrix.elementsEqual(mat(L(-2), L(1), L(4), L(-3)), Matrix.multiply(ssMat, ssTr)));

        Object[][] sq = Matrix.multiply(ssMat, ssMat);
        Object[][] expected = mat(L(1 + 6), L(-2 - 8), L(-3 - 12), L(6 + 16));
        assertTrue(Matrix.elementsEqual(expected, sq));
    }

    @Test
    public void testDiffSize() throws InitException {
        Object[][] m1 = mat(L(1), L(0), L(0), L(1));
        // Build a non-square matrix manually.
        Object[][] m2 = new Object[][] {
            new Object[] { L(1), L(0), L(0) },
            new Object[] { L(0), L(1) },
            new Object[] { L(0), L(0), L(1) },
        };
        try {
            Matrix.multiply(m1, m2);
            fail("Expected IllegalArgumentException for non-square");
        } catch (IllegalArgumentException ex) {
            // ok
        }
    }


    /* ----------------------- TestLinear ----------------------- */

    @Test
    public void testLinearIntSize2() throws InitException {
        Object[][] ssMat = mat(L(1), L(-2), L(-3), L(4));
        assertTrue(Matrix.elementsEqual(ssMat, Matrix.linear(ssMat)));
        assertTrue(Matrix.elementsEqual(mat(L(2), L(-4), L(-6), L(8)),
                                         Matrix.linear(ssMat, L(2))));
        assertTrue(Matrix.elementsEqual(mat(L(1), L(-5), L(-7), L(7)),
                                         Matrix.linear(ssMat, L(2), L(-1))));
    }

    @Test
    public void testLinearVardblSize2() throws InitException {
        Object[][] ssMat = new Object[][] {
            { new VarDbl(1.0, 1e-3), new VarDbl(2.0, 1e-2) },
            { new VarDbl(3.0), new VarDbl(4.0, 0.1) },
        };
        Object[][] r = Matrix.linear(ssMat, L(-2), L(1));
        verifyValue(new VarDbl(-1.0, 2e-3), r[0][0]);
        verifyValue(new VarDbl(-3.0, 2e-2), r[0][1]);
        verifyValue(new VarDbl(-5.0, 0), r[1][0]);
        verifyValue(new VarDbl(-7.0, 0.2), r[1][1]);
    }


    /* ----------------------- TestAdjugate ----------------------- */

    @Test
    public void testAdjIntSize2() throws InitException {
        verify(mat(L(1), L(0), L(0), L(1)), L(1), mat(L(1), L(0), L(0), L(1)));
        verify(mat(L(0), L(1), L(1), L(0)), L(-1), mat(L(0), L(-1), L(-1), L(0)));
        verify(mat(L(1), L(2), L(3), L(4)), L(-2), mat(L(4), L(-2), L(-3), L(1)));
    }

    @Test
    public void testAdjFloatSize2() throws InitException {
        Object[][] ssMat = mat(L(1), Double.valueOf(2.0), L(3), L(4));
        Object[][] ssAdj = mat(L(4), Double.valueOf(-2.0), L(-3), L(1));
        verify(ssMat, Double.valueOf(-2.0), ssAdj);
    }

    @Test
    public void testAdjIntFractionSize2() throws InitException {
        Fraction half = new Fraction(1, 2);
        Object[][] ssMat = mat(L(1), L(2), L(3), half);
        Object[][] ssAdj = mat(half, L(-2), L(-3), L(1));
        verify(ssMat, new Fraction(-11, 2), ssAdj);
    }

    @Test
    public void testAdjIntFractionFloatSize2() throws InitException {
        Fraction half = new Fraction(1, 2);
        Object[][] ssMat = mat(L(1), Double.valueOf(2.0), L(3), half);
        Object[][] ssAdj = mat(half, Double.valueOf(-2.0), L(-3), L(1));
        verify(ssMat, Double.valueOf(-11.0 / 2.0), ssAdj);
    }

    @Test
    public void testAdjVarDblSize2() throws InitException {
        Object[][] ssMat = new Object[][] {
            { new VarDbl(1.0, 1e-1), new VarDbl(2.0, 1e-2) },
            { new VarDbl(3.0, 1e-3), new VarDbl(4.0, 1e-4) },
        };
        VarDbl det = (VarDbl) new VarDbl(1.0, 1e-1).multiply(new VarDbl(4.0, 1e-4))
                .minus(new VarDbl(2.0, 1e-2).multiply(new VarDbl(3.0, 1e-3)));
        Object[][] ssAdj = new Object[][] {
            { new VarDbl(4.0, 1e-4), new VarDbl(-2.0, 1e-2) },
            { new VarDbl(-3.0, 1e-3), new VarDbl(1.0, 1e-1) },
        };
        verify(ssMat, det, ssAdj);
    }

    @Test
    public void testAdjIntSize3() throws InitException {
        verify(mat(L(1), L(0), L(0), L(0), L(1), L(0), L(0), L(0), L(1)),
               L(1),
               mat(L(1), L(0), L(0), L(0), L(1), L(0), L(0), L(0), L(1)));
        verify(mat(L(1), L(2), L(3), L(-4), L(-5), L(6), L(7), L(8), L(9)),
               L(72),
               mat(L(-93), L(6), L(27), L(78), L(-12), L(-18), L(3), L(6), L(3)));
    }

    @Test
    public void testVarDblSize3_noTracing() throws InitException {
        // Same variance computation as the Python test.
        VarDbl noTracing =
            (VarDbl)((VarDbl) new VarDbl(1.0, 1e-1).multiply(new VarDbl(-5.0, 1e-5)).multiply(new VarDbl(9.0, 1e-9)))
            .minus((VarDbl) new VarDbl(1.0, 1e-1).multiply(new VarDbl(6.0, 1e-6)).multiply(new VarDbl(8.0, 1e-8)))
            .add((VarDbl) new VarDbl(2.0, 1e-2).multiply(new VarDbl(6.0, 1e-6)).multiply(new VarDbl(7.0, 1e-7)))
            .minus((VarDbl) new VarDbl(2.0, 1e-2).multiply(new VarDbl(-4.0, 1e-4)).multiply(new VarDbl(9.0, 1e-9)))
            .add((VarDbl) new VarDbl(3.0, 1e-3).multiply(new VarDbl(-4.0, 1e-4)).multiply(new VarDbl(8.0, 1e-8)))
            .minus((VarDbl) new VarDbl(3.0, 1e-3).multiply(new VarDbl(-5.0, 1e-5)).multiply(new VarDbl(7.0, 1e-7)));
        verifyValue(new VarDbl(72, 6.602897701208691), noTracing);
    }

    @Test
    public void testAdjVarDblSize3() throws InitException {
        Object[][] ssMat3 = new Object[][] {
            { new VarDbl(1.0, 1e-1), new VarDbl(2.0, 1e-2), new VarDbl(3.0, 1e-3) },
            { new VarDbl(-4.0, 1e-4), new VarDbl(-5.0, 1e-5), new VarDbl(6.0, 1e-6) },
            { new VarDbl(7.0, 1e-7), new VarDbl(8.0, 1e-8), new VarDbl(9.0, 1e-9) },
        };
        double val3 = 72;
        double var3 =
              Math.pow(1e-1 * 1e-5 * 1e-9, 2) + Math.pow(1e-1 * 1e-6 * 1e-8, 2)
            + Math.pow(1e-2 * 1e-6 * 1e-7, 2) + Math.pow(1e-2 * 1e-4 * 1e-9, 2)
            + Math.pow(1e-3 * 1e-4 * 1e-8, 2) + (1e-3 * 1e-5 * 1e-7)
            + (1 * 1) * 1e-5 * 1e-9 + 1e-1 * (5 * 5) * 1e-9 + 1e-1 * 1e-5 * (9 * 9)
            + (1 * 1) * 1e-6 * 1e-8 + 1e-1 * (6 * 6) * 1e-8 + 1e-1 * 1e-6 * (8 * 8)
            + (2 * 2) * 1e-6 * 1e-7 + 1e-2 * (6 * 6) * 1e-7 + 1e-2 * 1e-6 * (7 * 7)
            + (2 * 2) * 1e-4 * 1e-9 + 1e-2 * (4 * 4) * 1e-9 + 1e-2 * 1e-4 * (9 * 9)
            + (3 * 3) * 1e-4 * 1e-8 + 1e-3 * (4 * 4) * 1e-8 + 1e-3 * 1e-4 * (8 * 8)
            + (3 * 3) * 1e-5 * 1e-7 + 1e-3 * (5 * 5) * 1e-7 + 1e-3 * 1e-5 * (7 * 7)
            + Math.pow(1e-1 * ((-5) * 9 - 6 * 8), 2)
            + Math.pow(1e-2 * ((-4) * 9 - 6 * 7), 2)
            + Math.pow(1e-3 * ((-4) * 8 - (-5) * 7), 2)
            + Math.pow(1e-4 * (2 * 9 - 3 * 8), 2)
            + Math.pow(1e-5 * (1 * 9 - 3 * 7), 2)
            + Math.pow(1e-6 * (1 * 8 - 2 * 7), 2)
            + Math.pow(1e-7 * (2 * 6 - 3 * (-5)), 2)
            + Math.pow(1e-8 * (1 * 6 - 3 * (-4)), 2)
            + Math.pow(1e-9 * (1 * (-5) - 2 * (-4)), 2);
        assertEquals(87.09858523178215, var3, 1e-9);

        Object[][] ssAdj = new Object[][] {
            {
                new VarDbl(-93, Math.sqrt(Math.pow(1e-5 * 1e-9, 2) + Math.pow(5 * 1e-9, 2) + Math.pow(1e-5 * 9, 2)
                                          + Math.pow(1e-6 * 1e-8, 2) + Math.pow(6 * 1e-8, 2) + Math.pow(1e-6 * 8, 2))),
                new VarDbl(6, Math.sqrt(Math.pow(1e-2 * 1e-9, 2) + Math.pow(2 * 1e-9, 2) + Math.pow(1e-2 * 9, 2)
                                          + Math.pow(1e-8 * 1e-3, 2) + Math.pow(8 * 1e-3, 2) + Math.pow(1e-8 * 3, 2))),
                new VarDbl(27, Math.sqrt(Math.pow(1e-2 * 1e-6, 2) + Math.pow(2 * 1e-5, 2) + Math.pow(1e-2 * 6, 2)
                                          + Math.pow(1e-5 * 1e-3, 2) + Math.pow(5 * 1e-3, 2) + Math.pow(1e-5 * 3, 2))),
            },
            {
                new VarDbl(78, Math.sqrt(Math.pow(1e-5 * 1e-7, 2) + Math.pow(6 * 1e-7, 2) + Math.pow(1e-5 * 7, 2)
                                          + Math.pow(1e-4 * 1e-9, 2) + Math.pow(4 * 1e-9, 2) + Math.pow(1e-4 * 9, 2))),
                new VarDbl(-12, Math.sqrt(Math.pow(1e-1 * 1e-9, 2) + Math.pow(1 * 1e-9, 2) + Math.pow(1e-1 * 9, 2)
                                          + Math.pow(1e-3 * 1e-7, 2) + Math.pow(3 * 1e-7, 2) + Math.pow(1e-3 * 7, 2))),
                new VarDbl(-18, Math.sqrt(Math.pow(1e-1 * 1e-6, 2) + Math.pow(1 * 1e-6, 2) + Math.pow(1e-1 * 6, 2)
                                          + Math.pow(1e-4 * 1e-3, 2) + Math.pow(4 * 1e-3, 2) + Math.pow(1e-4 * 3, 2))),
            },
            {
                new VarDbl(3, Math.sqrt(Math.pow(1e-4 * 1e-8, 2) + Math.pow(4 * 1e-8, 2) + Math.pow(1e-4 * 8, 2)
                                          + Math.pow(1e-5 * 1e-7, 2) + Math.pow(5 * 1e-7, 2) + Math.pow(1e-5 * 7, 2))),
                new VarDbl(6, Math.sqrt(Math.pow(1e-2 * 1e-7, 2) + Math.pow(2 * 1e-7, 2) + Math.pow(1e-2 * 7, 2)
                                          + Math.pow(1e-1 * 1e-8, 2) + Math.pow(1 * 1e-8, 2) + Math.pow(1e-1 * 8, 2))),
                new VarDbl(3, Math.sqrt(Math.pow(1e-1 * 1e-5, 2) + Math.pow(1 * 1e-5, 2) + Math.pow(1e-1 * 5, 2)
                                          + Math.pow(1e-2 * 1e-4, 2) + Math.pow(2 * 1e-4, 2) + Math.pow(1e-2 * 4, 2))),
            },
        };
        verify(ssMat3, new VarDbl(val3, Math.sqrt(var3)), ssAdj);
    }

    @Test
    public void testDirectMultiplication() throws InitException {
        Object[][] ssMat3 = new Object[][] {
            { new VarDbl(1.0, 1e-1), new VarDbl(2.0, 1e-2), new VarDbl(3.0, 1e-3) },
            { new VarDbl(-4.0, 1e-4), new VarDbl(-5.0, 1e-5), new VarDbl(6.0, 1e-6) },
            { new VarDbl(7.0, 1e-7), new VarDbl(8.0, 1e-8), new VarDbl(9.0, 1e-9) },
        };
        Matrix.AdjugateResult r = Matrix.adjugate_mul(ssMat3);
        assertTrue("det should be VarDbl, got " + r.det.getClass(), r.det instanceof VarDbl);
        assertEquals(72.0, ((VarDbl) r.det).value(), 1e-7);
        long[] expected = { -93, 6, 27, 78, -12, -18, 3, 6, 3 };
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Object cell = r.adj[i][j];
                double v = (cell instanceof VarDbl) ? ((VarDbl) cell).value()
                            : ((Number) cell).doubleValue();
                assertEquals(String.format("adj[%d][%d]", i, j), expected[i * 3 + j], v, 1e-7);
            }
        }
    }


    /* ----------------------- testInversionValue ----------------------- */

    @Test
    public void testInversionValue() throws InitException {
        // Pure scalar arithmetic check from Python.
        double common =
            ((1 * 1) * Math.pow(1e-4, 2) + (4 * 4) * Math.pow(1e-1, 2)
             + (2 * 2) * Math.pow(1e-3, 2) + (3 * 3) * Math.pow(1e-2, 2))
            / Math.pow(2, 4);
        assertEquals(0.010056500625000003, common, 1e-15);

        double[][] sX = {
            { 4.0 / -2 - 1.0 / (2 * 2) * Math.pow(1e-4, 2),
              -2.0 / -2 - -3.0 / (2 * 2) * Math.pow(1e-2, 2) },
            { -3.0 / -2 - -2.0 / (2 * 2) * Math.pow(1e-3, 2),
              1.0 / -2 - 4.0 / (2 * 2) * Math.pow(1e-1, 2) },
        };
        assertEquals(-2.0000000025, sX[0][0], 1e-12);
        assertEquals(1.000075, sX[0][1], 1e-12);
        assertEquals(1.5000005, sX[1][0], 1e-12);
        assertEquals(-0.51, sX[1][1], 1e-12);

        double[][] sdX = {
            { Math.sqrt(Math.pow(1e-4, 2) * ((1.0 / (2 * 2)) + 4.0 * 1 * 4 / Math.pow(2, 3))
                        + (4 * 4) * common),
              Math.sqrt(Math.pow(1e-2, 2) * ((1.0 / (2 * 2)) + 4.0 * 2 * 3 / Math.pow(2, 3))
                        + (2 * 2) * common) },
            { Math.sqrt(Math.pow(1e-3, 2) * ((1.0 / (2 * 2)) + 4.0 * 2 * 3 / Math.pow(2, 3))
                        + (3 * 3) * common),
              Math.sqrt(Math.pow(1e-1, 2) * ((1.0 / (2 * 2)) + 4.0 * 1 * 4 / Math.pow(2, 3))
                        + (1 * 1) * common) },
        };
        assertEquals(0.40112844887890964, sdX[0][0], 1e-15);
        assertEquals(0.2013727948358467, sdX[0][1], 1e-15);
        assertEquals(0.30085171700523833, sdX[1][0], 1e-15);
        assertEquals(0.1804342002642515, sdX[1][1], 1e-15);
    }
}
