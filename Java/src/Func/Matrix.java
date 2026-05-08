package Func;

/**
 * Matrix utilities for VarDbl/Fraction/Double/Long element types: permutation
 * sign, square-matrix checks, integer/Hilbert matrix generation, linear scale,
 * multiply, adjugate, and adjugate_mul. Mirrors Python matrix.py.
 *
 * Element types are carried as boxed Object instances so a single matrix may
 * hold a mix of integers (Long), rationals (Fraction), doubles (Double), and
 * VarDbl. The promotion order in mixed arithmetic is Long -> Fraction ->
 * Double -> VarDbl, matching the Python implementation.
 */
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.fraction.Fraction;

import Type.InitException;
import Type.VarDbl;


public class Matrix {

    /**
     * Tuple-style return of (det, adj) from {@link #adjugate(Object[][])}.
     */
    public static class AdjugateResult {
        public final Object det;
        public final Object[][] adj;
        public AdjugateResult(Object det, Object[][] adj) {
            this.det = det;
            this.adj = adj;
        }
    }


    /**
     * +1 or -1 according to the parity of the permutation.
     * https://en.wikipedia.org/wiki/Parity_of_a_permutation
     */
    public static int permutSign(int[] sPermut) {
        int cnt = 0;
        for (int i = 0; i < sPermut.length - 1; ++i) {
            for (int j = i + 1; j < sPermut.length; ++j) {
                if (sPermut[i] > sPermut[j]) {
                    cnt += 1;
                }
            }
        }
        return ((cnt % 2) != 0) ? -1 : +1;
    }


    /**
     * Generate every permutation of {0,...,size-1}. Lexicographic order
     * matching Python's itertools.permutations.
     */
    public static List<int[]> permutations(int size) {
        List<int[]> result = new ArrayList<int[]>();
        int[] state = new int[size];
        for (int i = 0; i < size; ++i) state[i] = i;
        boolean[] used = new boolean[size];
        int[] current = new int[size];
        permute(0, size, used, current, result);
        return result;
    }

    private static void permute(int depth, int size, boolean[] used, int[] current, List<int[]> out) {
        if (depth == size) {
            out.add(current.clone());
            return;
        }
        for (int i = 0; i < size; ++i) {
            if (used[i]) continue;
            used[i] = true;
            current[depth] = i;
            permute(depth + 1, size, used, current, out);
            used[i] = false;
        }
    }


    /**
     * Generate every k-combination of {0,...,n-1}. Lexicographic order
     * matching Python's itertools.combinations.
     */
    public static List<int[]> combinations(int n, int k) {
        List<int[]> result = new ArrayList<int[]>();
        if ((k < 0) || (k > n)) return result;
        int[] cur = new int[k];
        combine(0, 0, n, k, cur, result);
        return result;
    }

    private static void combine(int start, int depth, int n, int k, int[] cur, List<int[]> out) {
        if (depth == k) {
            out.add(cur.clone());
            return;
        }
        for (int i = start; i < n; ++i) {
            cur[depth] = i;
            combine(i + 1, depth + 1, n, k, cur, out);
        }
    }


    private static boolean isAllowedType(Object v) {
        return (v instanceof Long) || (v instanceof Integer)
            || (v instanceof Double) || (v instanceof Float)
            || (v instanceof Fraction) || (v instanceof VarDbl);
    }

    /**
     * Match Python isSquareMatrix. Java arrays are always rectangular, but we
     * still verify length and element types (rejecting nulls / unknown types).
     * If "allowedTypes" is non-null, every element must match one of the listed
     * classes.
     */
    public static boolean isSquareMatrix(Object ssMatrix, Class<?>[] allowedTypes) {
        if (!(ssMatrix instanceof Object[][])) return false;
        Object[][] mat = (Object[][]) ssMatrix;
        int size = mat.length;
        for (Object[] row : mat) {
            if (row == null || row.length != size) return false;
            for (Object v : row) {
                if (v == null) return false;
                if (allowedTypes != null) {
                    boolean ok = false;
                    for (Class<?> c : allowedTypes) {
                        if (c.isInstance(v)) { ok = true; break; }
                    }
                    if (!ok) return false;
                } else if (!isAllowedType(v)) {
                    return false;
                }
            }
        }
        return true;
    }

    public static boolean isSquareMatrix(Object ssMatrix) {
        return isSquareMatrix(ssMatrix, null);
    }


    /**
     * Random int matrix of the requested size with elements drawn uniformly
     * from [-randRange, +randRange].
     */
    public static Object[][] createIntMatrix(int size, int randRange, java.util.Random rng) {
        size = Math.abs(size);
        Object[][] out = new Object[size][size];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                int v = rng.nextInt(2 * randRange + 1) - randRange;
                out[i][j] = Long.valueOf(v);
            }
        }
        return out;
    }


    /**
     * Hilbert matrix as Fractions. https://en.wikipedia.org/wiki/Hilbert_matrix
     */
    public static Object[][] createHilbertMatrix(int size) {
        size = Math.abs(size);
        Object[][] out = new Object[size][size];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                out[i][j] = new Fraction(1, i + j + 1);
            }
        }
        return out;
    }


    /* ===================== Element-wise arithmetic helpers ===================== */

    private static double toDouble(Object v) {
        if (v instanceof Long) return ((Long) v).doubleValue();
        if (v instanceof Integer) return ((Integer) v).doubleValue();
        if (v instanceof Fraction) return ((Fraction) v).doubleValue();
        if (v instanceof Double) return ((Double) v).doubleValue();
        if (v instanceof Float) return ((Float) v).doubleValue();
        if (v instanceof VarDbl) return ((VarDbl) v).value();
        throw new IllegalArgumentException("Unsupported element type: " + v);
    }

    /** Promotion rank: Long < Fraction < Double < VarDbl. */
    private static int rank(Object v) {
        if (v instanceof Long || v instanceof Integer) return 0;
        if (v instanceof Fraction) return 1;
        if (v instanceof Double || v instanceof Float) return 2;
        if (v instanceof VarDbl) return 3;
        throw new IllegalArgumentException("Unsupported element type: " + v);
    }

    public static Object negate(Object a) {
        if (a instanceof Long) return Long.valueOf(-((Long) a).longValue());
        if (a instanceof Integer) return Long.valueOf(-((Integer) a).longValue());
        if (a instanceof Fraction) return ((Fraction) a).negate();
        if (a instanceof Double) return Double.valueOf(-((Double) a).doubleValue());
        if (a instanceof Float) return Double.valueOf(-((Float) a).doubleValue());
        if (a instanceof VarDbl) {
            try {
                return new VarDbl((VarDbl) a).negate();
            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
        throw new IllegalArgumentException("Unsupported element type: " + a);
    }

    /** Zero of the dominant rank type. */
    private static Object zeroOfRank(int rk) {
        switch (rk) {
            case 0: return Long.valueOf(0L);
            case 1: return Fraction.ZERO;
            case 2: return Double.valueOf(0.0);
            case 3:
                try { return new VarDbl(0.0, 0.0); }
                catch (InitException ex) { throw new RuntimeException(ex); }
            default:
                throw new IllegalArgumentException("rank " + rk);
        }
    }

    public static Object add(Object a, Object b) throws InitException {
        int ra = rank(a), rb = rank(b);
        int r = Math.max(ra, rb);
        if (r == 0) {
            long la = (a instanceof Long) ? (Long) a : ((Integer) a).longValue();
            long lb = (b instanceof Long) ? (Long) b : ((Integer) b).longValue();
            return Long.valueOf(la + lb);
        }
        if (r == 1) {
            Fraction fa = toFraction(a);
            Fraction fb = toFraction(b);
            return fa.add(fb);
        }
        if (r == 2) {
            return Double.valueOf(toDouble(a) + toDouble(b));
        }
        // VarDbl
        VarDbl va = toVarDbl(a);
        VarDbl vb = toVarDbl(b);
        return va.add(vb);
    }

    public static Object multiply(Object a, Object b) throws InitException {
        int ra = rank(a), rb = rank(b);
        int r = Math.max(ra, rb);
        if (r == 0) {
            long la = (a instanceof Long) ? (Long) a : ((Integer) a).longValue();
            long lb = (b instanceof Long) ? (Long) b : ((Integer) b).longValue();
            return Long.valueOf(la * lb);
        }
        if (r == 1) {
            Fraction fa = toFraction(a);
            Fraction fb = toFraction(b);
            return fa.multiply(fb);
        }
        if (r == 2) {
            return Double.valueOf(toDouble(a) * toDouble(b));
        }
        VarDbl va = toVarDbl(a);
        VarDbl vb = toVarDbl(b);
        return va.multiply(vb);
    }

    private static Fraction toFraction(Object v) {
        if (v instanceof Long) {
            long lv = (Long) v;
            if (lv > Integer.MAX_VALUE || lv < Integer.MIN_VALUE)
                throw new ArithmeticException("Long out of int range for Fraction: " + lv);
            return new Fraction((int) lv, 1);
        }
        if (v instanceof Integer) return new Fraction((Integer) v, 1);
        if (v instanceof Fraction) return (Fraction) v;
        throw new IllegalArgumentException("Cannot convert to Fraction: " + v);
    }

    private static VarDbl toVarDbl(Object v) throws InitException {
        if (v instanceof VarDbl) return (VarDbl) v;
        if (v instanceof Long) return new VarDbl(((Long) v).longValue());
        if (v instanceof Integer) return new VarDbl(((Integer) v).longValue());
        if (v instanceof Fraction) return new VarDbl(((Fraction) v).doubleValue(), 0.0);
        if (v instanceof Double) return new VarDbl(((Double) v).doubleValue());
        if (v instanceof Float) return new VarDbl(((Float) v).floatValue());
        throw new IllegalArgumentException("Cannot convert to VarDbl: " + v);
    }


    /* ============================ Linear ============================ */

    /** linear(M, scale, offset) = M*scale + offset element-wise. */
    public static Object[][] linear(Object[][] ssMatrix, Object scale, Object offset) throws InitException {
        if (!isSquareMatrix(ssMatrix))
            throw new IllegalArgumentException("Invalid square matrix for linear()");
        if (!isAllowedType(scale))
            throw new IllegalArgumentException("Invalid scale type for linear(): " + scale);
        int size = ssMatrix.length;
        Object[][] out = new Object[size][size];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Object scaled = multiply(ssMatrix[i][j], scale);
                Object res = add(scaled, offset);
                out[i][j] = res;
            }
        }
        return out;
    }

    public static Object[][] linear(Object[][] ssMatrix, Object scale) throws InitException {
        return linear(ssMatrix, scale, Long.valueOf(0L));
    }

    public static Object[][] linear(Object[][] ssMatrix) throws InitException {
        return linear(ssMatrix, Long.valueOf(1L), Long.valueOf(0L));
    }


    /* ============================ Multiply ============================ */

    public static Object[][] multiply(Object[][] a, Object[][] b) throws InitException {
        if (!isSquareMatrix(a))
            throw new IllegalArgumentException("Invalid square matrix 1 for multiply()");
        if (!isSquareMatrix(b))
            throw new IllegalArgumentException("Invalid square matrix 2 for multiply()");
        int size = a.length;
        if (size != b.length)
            throw new IllegalArgumentException(
                String.format("Different sizes for multiply(): %d vs %d", size, b.length));
        Object[][] out = new Object[size][size];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                // Determine rank by inspection so we sum into the correct accumulator.
                int maxR = 0;
                for (int k = 0; k < size; ++k) {
                    maxR = Math.max(maxR, Math.max(rank(a[i][k]), rank(b[k][j])));
                }
                Object acc = zeroOfRank(maxR);
                for (int k = 0; k < size; ++k) {
                    Object prod = multiply(a[i][k], b[k][j]);
                    acc = add(acc, prod);
                }
                out[i][j] = acc;
            }
        }
        return out;
    }


    /* ============================ Adjugate ============================ */

    /**
     * Determinant and adjugate. Mirrors Python adjugate(). Promotion order is
     * Long -> Fraction -> Double -> VarDbl. For VarDbl input the variance is
     * propagated through the variance-arithmetic algorithm in matrix.py.
     */
    public static AdjugateResult adjugate(Object[][] ssMatrix) throws InitException {
        if (!isSquareMatrix(ssMatrix))
            throw new IllegalArgumentException("Invalid square matrix for adjugate()");
        int size = ssMatrix.length;
        if (size == 1) {
            Object[][] adj = new Object[1][1];
            adj[0][0] = ssMatrix[0][0];
            return new AdjugateResult(ssMatrix[0][0], adj);
        }
        if (size == 2) {
            Object det = add(multiply(ssMatrix[0][0], ssMatrix[1][1]),
                            negate(multiply(ssMatrix[0][1], ssMatrix[1][0])));
            Object[][] adj = new Object[2][2];
            adj[0][0] = ssMatrix[1][1];
            adj[0][1] = negate(ssMatrix[0][1]);
            adj[1][0] = negate(ssMatrix[1][0]);
            adj[1][1] = ssMatrix[0][0];
            return new AdjugateResult(det, adj);
        }

        // size >= 3 — use the variance-tracing algorithm.
        // We carry value/variance separately. Determine top rank: VarDbl if any,
        // else fall back to elementwise product types.
        boolean anyVarDbl = false;
        int topRank = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (ssMatrix[i][j] instanceof VarDbl) anyVarDbl = true;
                topRank = Math.max(topRank, rank(ssMatrix[i][j]));
            }
        }

        List<int[]> sPermut = permutations(size);
        int[] sSign = new int[sPermut.size()];
        for (int p = 0; p < sPermut.size(); ++p) {
            sSign[p] = permutSign(sPermut.get(p));
        }

        // Accumulator for determinant value and variance.
        Object value = zeroOfRank(topRank);
        double variance = 0.0;
        // Cofactor variance map keyed by (x,y) packed as x*size+y.
        double[] sCofVar = new double[size * size];

        for (int p = 0; p < sPermut.size(); ++p) {
            int[] permut = sPermut.get(p);
            int sign = sSign[p];
            // val starts as scaled by sign (long) and grows in rank.
            Object val = Long.valueOf(sign);
            double var = 1.0;
            boolean hasVarDbl = false;
            for (int x = 0; x < size; ++x) {
                int y = permut[x];
                Object e = ssMatrix[x][y];
                if (e instanceof VarDbl) {
                    hasVarDbl = true;
                    val = multiply(val, Double.valueOf(((VarDbl) e).value()));
                    var *= ((VarDbl) e).variance();
                } else {
                    val = multiply(val, e);
                    var = 0.0;
                }
            }
            if (!hasVarDbl) var = 0.0;
            value = add(value, val);
            variance += var;
            if (var > 0) {
                for (int x = 0; x < size; ++x) {
                    int y = permut[x];
                    if (ssMatrix[x][y] instanceof VarDbl) {
                        double v = ((VarDbl) ssMatrix[x][y]).variance();
                        if (v > 0) {
                            sCofVar[x * size + y] += var / v;
                        }
                    }
                }
            }
        }

        // Per-cofactor value tracking. sCofVal[x*size+y] holds the (signed) cofactor value.
        // We initialise lazily (m == 1 step builds it).
        Object[] sCofVal = new Object[size * size];

        for (int m = 1; m < size; ++m) {
            // Map keyed by sY tuple-of-pairs (encoded as a string).
            HashMap<String, Object> sVal = new HashMap<String, Object>();
            HashMap<String, Double> sVar = new HashMap<String, Double>();
            HashMap<String, int[][]> sYMap = new HashMap<String, int[][]>();

            List<int[]> combos = combinations(size, m);
            for (int[] sX : combos) {
                for (int p = 0; p < sPermut.size(); ++p) {
                    int[] permut = sPermut.get(p);
                    int sign = sSign[p];
                    Object val = Long.valueOf(sign);
                    double var = 1.0;
                    int[][] sY = new int[m][2];
                    for (int xi = 0; xi < m; ++xi) {
                        int xx = sX[xi];
                        sY[xi][0] = xx;
                        sY[xi][1] = permut[xx];
                    }
                    // Walk the permutation.
                    int xi = 0;
                    for (int x = 0; x < size; ++x) {
                        int y = permut[x];
                        boolean inSX = (xi < m) && (sX[xi] == x);
                        if (inSX) {
                            ++xi;
                            if (ssMatrix[x][y] instanceof VarDbl) {
                                var *= ((VarDbl) ssMatrix[x][y]).variance();
                            } else {
                                var *= 0.0;
                            }
                        } else {
                            if (ssMatrix[x][y] instanceof VarDbl) {
                                val = multiply(val, Double.valueOf(((VarDbl) ssMatrix[x][y]).value()));
                            } else {
                                val = multiply(val, ssMatrix[x][y]);
                            }
                        }
                    }
                    String key = encodeSY(sY);
                    if (sVal.containsKey(key)) {
                        sVal.put(key, add(sVal.get(key), val));
                        // assert sVar[key] == var (by construction it should be — same sX has same product of variances)
                    } else {
                        sVal.put(key, val);
                        sVar.put(key, Double.valueOf(var));
                        sYMap.put(key, sY);
                    }
                }
            }

            for (Map.Entry<String, Double> entry : sVar.entrySet()) {
                String key = entry.getKey();
                double var = entry.getValue();
                Object vv = sVal.get(key);
                double vd = toDouble(vv);
                double inc = var * vd * vd;
                variance += inc;
                if (var > 0) {
                    int[][] sY = sYMap.get(key);
                    for (int[] xy : sY) {
                        int x = xy[0], y = xy[1];
                        if (ssMatrix[x][y] instanceof VarDbl) {
                            double evar = ((VarDbl) ssMatrix[x][y]).variance();
                            if (evar > 0) {
                                sCofVar[x * size + y] += var * vd * vd / evar;
                            }
                        }
                    }
                }
            }
            if (m == 1) {
                // sVal keys hold a single (x,y) pair; populate sCofVal[x*size+y].
                for (Map.Entry<String, Object> e : sVal.entrySet()) {
                    int[][] sY = sYMap.get(e.getKey());
                    int x = sY[0][0], y = sY[0][1];
                    sCofVal[x * size + y] = e.getValue();
                    double vd = toDouble(e.getValue());
                    sCofVar[x * size + y] -= vd * vd;
                }
            }
        }

        // Build result determinant.
        Object det;
        if (anyVarDbl && variance > 0) {
            double dv = toDouble(value);
            det = new VarDbl(dv, Math.sqrt(variance));
        } else {
            det = value;
        }

        // Build result adjugate matrix: out[j][i] = sCofVal[(i,j)] (note transpose)
        Object[][] adj = new Object[size][size];
        for (int j = 0; j < size; ++j) {
            for (int i = 0; i < size; ++i) {
                Object cv = sCofVal[i * size + j];
                double cVar = sCofVar[i * size + j];
                if (anyVarDbl && (cVar > 0)) {
                    double v = toDouble(cv);
                    adj[j][i] = new VarDbl(v, Math.sqrt(cVar));
                } else {
                    adj[j][i] = cv;
                }
            }
        }

        return new AdjugateResult(det, adj);
    }

    private static String encodeSY(int[][] sY) {
        StringBuilder sb = new StringBuilder();
        for (int[] xy : sY) {
            sb.append(xy[0]).append(',').append(xy[1]).append(';');
        }
        return sb.toString();
    }


    /* ============================ adjugate_mul ============================ */

    /**
     * Same result as {@link #adjugate(Object[][])} but computed by direct
     * multiplication (no variance tracing). Used to cross-check.
     */
    public static AdjugateResult adjugate_mul(Object[][] ssMatrix) throws InitException {
        if (!isSquareMatrix(ssMatrix))
            throw new IllegalArgumentException("Invalid square matrix for adjugate_mul()");
        int size = ssMatrix.length;

        // Determine zero rank.
        int topRank = 0;
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                topRank = Math.max(topRank, rank(ssMatrix[i][j]));

        Object det = zeroOfRank(topRank);
        Object[][] cof = new Object[size][size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cof[i][j] = zeroOfRank(topRank);

        List<int[]> sPermut = permutations(size);
        for (int[] permut : sPermut) {
            int sign = permutSign(permut);
            // Full product of entries along permutation.
            Object full = Long.valueOf(sign);
            for (int x = 0; x < size; ++x) {
                full = multiply(full, ssMatrix[x][permut[x]]);
            }
            det = add(det, full);

            // For each (i, permut[i]) the cofactor entry gets "sign * product over x != i".
            for (int i = 0; i < size; ++i) {
                Object prod = Long.valueOf(sign);
                for (int x = 0; x < size; ++x) {
                    if (x == i) continue;
                    prod = multiply(prod, ssMatrix[x][permut[x]]);
                }
                cof[i][permut[i]] = add(cof[i][permut[i]], prod);
            }
        }

        // Python returns adj[(i,j)] indexed as out[j][i] (transpose) just like adjugate.
        Object[][] adj = new Object[size][size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                adj[j][i] = cof[i][j];

        return new AdjugateResult(det, adj);
    }


    /* ============================ Public type helpers ============================ */

    /**
     * Convenience: convert a 2D array of doubles to Object[][] of Doubles.
     */
    public static Object[][] box(double[][] m) {
        int n = m.length;
        Object[][] out = new Object[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                out[i][j] = Double.valueOf(m[i][j]);
        return out;
    }

    /**
     * Convenience: convert a 2D array of longs to Object[][] of Longs.
     */
    public static Object[][] box(long[][] m) {
        int n = m.length;
        Object[][] out = new Object[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                out[i][j] = Long.valueOf(m[i][j]);
        return out;
    }

    public static Object[][] box(int[][] m) {
        int n = m.length;
        Object[][] out = new Object[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                out[i][j] = Long.valueOf(m[i][j]);
        return out;
    }

    /** Pretty-print helper, primarily for debugging tests. */
    public static String pretty(Object[][] m) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < m.length; ++i) {
            if (i > 0) sb.append(", ");
            sb.append(Arrays.toString(m[i]));
        }
        sb.append("]");
        return sb.toString();
    }

    /** True if a Object[][] m equals another using element-by-element semantic equality. */
    public static boolean elementsEqual(Object[][] a, Object[][] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; ++i) {
            if (a[i].length != b[i].length) return false;
            for (int j = 0; j < a[i].length; ++j) {
                if (!objectEquals(a[i][j], b[i][j])) return false;
            }
        }
        return true;
    }

    public static boolean objectEquals(Object a, Object b) {
        if (a == null || b == null) return a == b;
        if ((a instanceof Long || a instanceof Integer)
                && (b instanceof Long || b instanceof Integer)) {
            long la = (a instanceof Long) ? (Long) a : ((Integer) a).longValue();
            long lb = (b instanceof Long) ? (Long) b : ((Integer) b).longValue();
            return la == lb;
        }
        return a.equals(b);
    }
}
