package Func;

/**
 * Square matrix of VarDbl elements with determinant and adjugate computed via
 * bottom-up Laplace expansion for values (Formula 6.1) and direct evaluation
 * of Formula 6.9 for variances (independent-cell assumption).
 *
 * The static `runMatrixAnalysis` driver sweeps random + Hilbert matrices
 * across (size, noise) cells and writes MatrixCondition_*.txt and
 * AdjMatrix_*.txt under the given output directory. Cells are processed in
 * parallel via a fixed thread pool; the main thread drains results in
 * submission order so the output stays deterministic.
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import Stats.Stat;
import Type.InitException;
import Type.VarDbl;


public class Matrix {

    private final int size;
    private final VarDbl[] sValue;
    private final List<Set<Integer>> ssHigher = new ArrayList<>();

    /** Cached results — sAdjugate==null means "not yet computed". */
    private VarDbl determ;
    private VarDbl[] sAdjugate;

    /** Package-private — used internally to wrap an adjugate result. */
    Matrix(int size, VarDbl[] sValue) {
        this.size = size;
        this.sValue = sValue;
    }

    public Matrix(int size) throws InitException {
        this.size = size;
        this.sValue = new VarDbl[size * size];
        final VarDbl zero = new VarDbl();
        for (int i = 0; i < sValue.length; ++i)
             sValue[i] = zero;
    }

    public int    size() { return size; }

    private int index(int row, int col) {
        return row * size + col;
    }

    public VarDbl get(int row, int col) {
         return sValue[index(row, col)];
    }

    /** Set one cell; invalidates the cached adjugate/determinant. */
    public void set(int row, int col, VarDbl elem) {
        sValue[index(row, col)] = elem;
        sAdjugate = null;
    }
    /**
     * Set `elem` at every (row, col) in `positions` (must be non-empty),
     * invalidating any cached adjugate. If elem has positive variance AND the
     * caller supplied more than one position, the indices are recorded as a
     * higher-order correlated variance group; an index already in any prior
     * group fails the assertion.
     */
    public void set(VarDbl elem, List<int[]> positions) {
        assert positions != null && !positions.isEmpty()
            : "positions must be non-empty";
        sAdjugate = null;

        final Set<Integer> indices = new HashSet<>();
        final boolean grouped = elem.variance() > 0 && positions.size() > 1;
        for (int[] rc : positions) {
            assert rc != null && rc.length == 2 : "each position must be (row, col)";
            assert 0 <= rc[0] && rc[0] < size : "row out of range: " + rc[0];
            assert 0 <= rc[1] && rc[1] < size : "col out of range: " + rc[1];
            final int idx = index(rc[0], rc[1]);
            if (grouped) {
                for (Set<Integer> existing : ssHigher) {
                    assert !existing.contains(idx)
                        : "position " + idx + " already in a higher-order group";
                }
                indices.add(idx);
            }
            sValue[idx] = elem;
        }
        if (grouped) {
            ssHigher.add(indices);
        }
    }

    public VarDbl determ()    throws InitException {
         if (sAdjugate == null) calc();
          return determ;
    }

    public Matrix adjugate()  throws InitException {
        if (sAdjugate == null) calc();
        return new Matrix(size, sAdjugate.clone()); 
    }


    // =========================================================================
    //  calc(): determinant + adjugate via Formula 6.1 (value) and 6.9 (variance)
    // =========================================================================

    /**
     * Bottom-up Laplace builds a sub-determinant value table (doubles only),
     * then Formula 6.9 is evaluated directly for the determinant and for each
     * adjugate cell's submatrix:
     *
     *   δ²|M| = Σ_m Σ_{Sr, Sc both popcount m}
     *              |M_{~Sr,~Sc}|² · perm([σ²]_{Sr, Sc})
     *
     * The σ²-permanent collapses the ordered-tuple bijection sum into a single
     * recursive enumeration. Cost is Σ_m C(N,m)² · m! per submatrix — fine
     * through N ≈ 11, exponential beyond that.
     */
    private void calc() throws InitException {
        if (size == 0) {
            determ = new VarDbl();
            sAdjugate = new VarDbl[0];
            return;
        }
        final int nm = 1 << size;

        // Group masks by popcount; build mask→index lookup within each group.
        final int[][] masksAtCount = new int[size + 1][];
        {
            final int[] count = new int[size + 1];
            for (int m = 0; m < nm; ++m) count[Integer.bitCount(m)]++;
            for (int k = 0; k <= size; ++k) masksAtCount[k] = new int[count[k]];
            final int[] cursor = new int[size + 1];
            for (int m = 0; m < nm; ++m) {
                final int p = Integer.bitCount(m);
                masksAtCount[p][cursor[p]++] = m;
            }
        }
        final int[] maskIdx = new int[nm];
        for (int k = 0; k <= size; ++k)
            for (int i = 0; i < masksAtCount[k].length; ++i)
                maskIdx[masksAtCount[k][i]] = i;

        // Cell means and variances as primitive arrays.
        final double[] cellVal = new double[size * size];
        final double[] cellVar = new double[size * size];
        for (int idx = 0; idx < size * size; ++idx) {
            cellVal[idx] = sValue[idx].value();
            cellVar[idx] = sValue[idx].variance();
        }

        // Bottom-up Laplace: sub[k][rmIdx][cmIdx] = submatrix det at popcount k.
        final double[][][] sub = new double[size + 1][][];
        for (int k = 0; k <= size; ++k) {
            final int n = masksAtCount[k].length;
            sub[k] = new double[n][n];
        }
        sub[0][0][0] = 1.0;
        for (int k = 1; k <= size; ++k) {
            final int[] masks = masksAtCount[k];
            final double[][] subKm1 = sub[k - 1];
            for (int ai = 0; ai < masks.length; ++ai) {
                final int rm = masks[ai];
                final int firstRow = Integer.numberOfTrailingZeros(rm);
                final int rmRestIdx = maskIdx[rm & ~(1 << firstRow)];
                final int rowOff = firstRow * size;
                for (int bi = 0; bi < masks.length; ++bi) {
                    final int cm = masks[bi];
                    double det = 0.0;
                    int posJ = 0, rest = cm;
                    while (rest != 0) {
                        final int bit = rest & -rest;
                        final int j = Integer.numberOfTrailingZeros(bit);
                        final double term = cellVal[rowOff + j]
                                          * subKm1[rmRestIdx][maskIdx[cm ^ bit]];
                        det += ((posJ++ & 1) == 0) ? term : -term;
                        rest ^= bit;
                    }
                    sub[k][ai][bi] = det;
                }
            }
        }

        final int fullMask = nm - 1;
        final int fullIdx  = maskIdx[fullMask];

        // Determinant
        determ = new VarDbl(sub[size][fullIdx][fullIdx],
                            Math.sqrt(formula69(sub, masksAtCount, maskIdx,
                                                cellVar, fullMask, fullMask)));

        // Adjugate: cell value = ±(N−1)×(N−1) sub-det; variance from Formula 6.9.
        sAdjugate = new VarDbl[size * size];
        final double[][] subKm1Final = sub[size - 1];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                final int rmAdj = fullMask & ~(1 << j);
                final int cmAdj = fullMask & ~(1 << i);
                final double minor = subKm1Final[maskIdx[rmAdj]][maskIdx[cmAdj]];
                final double adjVal = ((i + j) & 1) == 0 ? minor : -minor;
                final double adjVar = formula69(sub, masksAtCount, maskIdx,
                                                cellVar, rmAdj, cmAdj);
                sAdjugate[i * size + j] = new VarDbl(adjVal, Math.sqrt(adjVar));
            }
        }
    }

    /** Formula 6.9 over the submatrix selected by (containerRowMask × containerColMask). */
    private double formula69(double[][][] sub, int[][] masksAtCount, int[] maskIdx,
                              double[] cellVar, int containerRowMask, int containerColMask) {
        final int nContainer = Integer.bitCount(containerRowMask);
        double variance = 0.0;
        for (int m = 1; m <= nContainer; ++m) {
            final int[] masks = masksAtCount[m];
            final int compLevel = nContainer - m;
            for (int posRowMask : masks) {
                if ((posRowMask & ~containerRowMask) != 0) continue;
                final int rmIdx = maskIdx[containerRowMask ^ posRowMask];
                for (int posColMask : masks) {
                    if ((posColMask & ~containerColMask) != 0) continue;
                    final double subVal = sub[compLevel][rmIdx][maskIdx[containerColMask ^ posColMask]];
                    final double sub2 = subVal * subVal;
                    if (sub2 == 0.0) continue;
                    final double perm = permanent(posRowMask, posColMask, cellVar);
                    if (perm > 0.0) variance += sub2 * perm;
                }
            }
        }
        return variance;
    }

    /** Permanent of the σ² submatrix selected by (rowsMask × colsMask). */
    private double permanent(int rowsMask, int colsMask, double[] cellVar) {
        if (rowsMask == 0) return 1.0;
        final int row = Integer.numberOfTrailingZeros(rowsMask);
        final int newRows = rowsMask & ~(1 << row);
        final int rowOff = row * size;
        double sum = 0.0;
        int rest = colsMask;
        while (rest != 0) {
            final int bit = rest & -rest;
            final int col = Integer.numberOfTrailingZeros(bit);
            rest ^= bit;
            final double sig = cellVar[rowOff + col];
            if (sig != 0.0) sum += sig * permanent(newRows, colsMask ^ bit, cellVar);
        }
        return sum;
    }


    // =========================================================================
    //  Multiplication
    // =========================================================================

    /** Returns this · other. Variance propagation under independent cells. */
    public Matrix multiply(final Matrix other) throws InitException {
        if (other.size != size)
            throw new IllegalArgumentException("size mismatch: " + size + " vs " + other.size);
        final VarDbl[] out = new VarDbl[size * size];
        for (int i = 0; i < size; ++i) {
            final int rowOff = i * size;
            for (int j = 0; j < size; ++j) {
                double accVal = 0.0, accVar = 0.0;
                for (int k = 0; k < size; ++k) {
                    final VarDbl a = sValue[rowOff + k];
                    final VarDbl b = other.sValue[k * size + j];
                    final double av = a.value(),       bv = b.value();
                    final double aVar = a.variance(), bVar = b.variance();
                    accVal += av * bv;
                    accVar += aVar * bv * bv + bVar * av * av + aVar * bVar;
                }
                out[i * size + j] = new VarDbl(accVal, Math.sqrt(accVar));
            }
        }
        return new Matrix(size, out);
    }


    // =========================================================================
    //  Static factories
    // =========================================================================

    /** N×N integer matrix with each cell uniform in [-range, range] inclusive. */
    public static int[][] createIntMatrix(int size, Random rng, int range) {
        final int[][] sInt = new int[size][size];
        final int span = 2 * range + 1;
        for (int r = 0; r < size; ++r)
            for (int c = 0; c < size; ++c) sInt[r][c] = rng.nextInt(span) - range;
        return sInt;
    }

    /** Wrap an integer matrix as a precise Matrix (cell variance = 0). */
    public static Matrix asMatrix(int[][] sInt) throws InitException {
        final int n = sInt.length;
        final VarDbl[] sv = new VarDbl[n * n];
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c) sv[r * n + c] = new VarDbl((long) sInt[r][c]);
        return new Matrix(n, sv);
    }

    /** Integer matrix with Gaussian noise of deviation `dev` added to each cell. */
    public static Matrix addNoise(int[][] sInt, double dev, Random rng) throws InitException {
        final int n = sInt.length;
        final VarDbl[] sv = new VarDbl[n * n];
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                sv[r * n + c] = new VarDbl(sInt[r][c] + dev * rng.nextGaussian(), dev);
        return new Matrix(n, sv);
    }

    /** Hilbert matrix H[i,j] = 1/(i+j+1), optionally noisy with deviation `dev`. */
    public static Matrix hilbertMatrix(int size, double dev, Random rng) throws InitException {
        final VarDbl[] sv = new VarDbl[size * size];
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                double v = 1.0 / (r + c + 1);
                if (dev > 0) v += dev * rng.nextGaussian();
                sv[r * size + c] = new VarDbl(v, dev);
            }
        }
        return new Matrix(size, sv);
    }

    /** Random matrix: each cell = (uniform int in [-128,128]) + Gaussian noise
     *  with deviation 256·dev; recorded uncertainty is 256·dev. */
    public static Matrix randomMatrix(int size, double dev) throws InitException {
        final Random rng = new Random();
        final double sigma = 256.0 * dev;
        final VarDbl[] sv = new VarDbl[size * size];
        for (int i = 0; i < sv.length; ++i) {
            final int u = rng.nextInt(257) - 128;
            sv[i] = new VarDbl(u + rng.nextGaussian() * sigma, sigma);
        }
        return new Matrix(size, sv);
    }


    // =========================================================================
    //  runMatrixAnalysis driver — port of Python/runMatrixAnalysis.py
    // =========================================================================

    private static final int    ELEMENT_RANGE = 1 << 8;
    private static final double DEV_SCALE     = ELEMENT_RANGE / Math.sqrt(3);

    /** Default 20-value noise sweep: {0, 1e-17, 1e-16, …, 1e1}. */
    public static double[] defaultNoises() {
        final double[] sN = new double[20];
        sN[0] = 0.0;
        for (int i = 1; i < 20; ++i) sN[i] = Math.pow(10, i - 18);
        return sN;
    }

    /** CLI args: [minSize] [maxSize] [targetCells] [outDir] [threads]. */
    public static void main(String[] args) throws InitException, IOException {
        final int minSize     = args.length > 0 ? Integer.parseInt(args[0]) : 4;
        final int maxSize     = args.length > 1 ? Integer.parseInt(args[1]) : 10;
        final int targetCells = args.length > 2 ? Integer.parseInt(args[2]) : 1000;
        final String outDir   = args.length > 3 ? args[3] : "Output";
        final int threads     = args.length > 4 ? Integer.parseInt(args[4])
                                                : Runtime.getRuntime().availableProcessors();
        System.out.printf("runMatrixAnalysis: sizes=[%d, %d) target=%d outDir=%s threads=%d%n",
                          minSize, maxSize, targetCells, outDir, threads);
        runMatrixAnalysis(minSize, maxSize, targetCells, defaultNoises(), outDir, threads);
        System.out.println("done");
    }

    /**
     * Sweep sizes × noises, accumulating ~targetCells cells per (size, noise).
     * Appends to existing output files and skips (size, noise) pairs already
     * present in the AdjMatrix file.
     */
    public static void runMatrixAnalysis(int minSize, int maxSize, int targetCells,
                                          double[] noises, String outDir, int threads)
            throws IOException {
        new File(outDir).mkdirs();
        final String condPath = outDir + "/MatrixCondition_" + minSize + "_" + maxSize + ".txt";
        final String adjPath  = outDir + "/AdjMatrix_"        + minSize + "_" + maxSize + ".txt";

        final Set<String> done = readDoneCells(adjPath);
        if (!done.isEmpty()) System.out.println("Skipping " + done.size() + " existing cells");

        final boolean condExists = new File(condPath).isFile();
        final boolean adjExists  = new File(adjPath).isFile();

        final ExecutorService exec = Executors.newFixedThreadPool(Math.max(1, threads));
        final List<Future<String[]>> futures = new ArrayList<>();
        try (BufferedWriter fc = new BufferedWriter(new FileWriter(condPath, condExists));
             BufferedWriter fa = new BufferedWriter(new FileWriter(adjPath,  adjExists))) {

            if (!condExists) fc.write(CONDITION_HEADER);
            if (!adjExists)  fa.write(ADJ_HEADER);

            for (int size = minSize; size < maxSize; ++size) {
                final int sz = size;
                final int count = Math.max(1, targetCells / (sz * sz));
                for (final double noise : noises) {
                    if (done.contains(sz + "|" + noise)) continue;
                    futures.add(exec.submit(() -> runCell(sz, noise, count)));
                }
            }
            exec.shutdown();

            for (Future<String[]> f : futures) {
                try {
                    final String[] r = f.get();
                    fc.write(r[0]); fc.flush();
                    fa.write(r[1]); fa.flush();
                } catch (ExecutionException ex) {
                    System.err.println("cell failed: " + ex.getCause());
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    throw new IOException("interrupted", ex);
                }
            }
        } finally {
            exec.shutdownNow();
        }
    }

    /** Read existing AdjMatrix file to determine which (size, noise) pairs are done. */
    private static Set<String> readDoneCells(String adjPath) throws IOException {
        final Set<String> done = new HashSet<>();
        final File f = new File(adjPath);
        if (!f.isFile()) return done;
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            br.readLine();      // skip header
            for (String line; (line = br.readLine()) != null; ) {
                final String[] parts = line.split("\t", -1);
                if (parts.length < 4) continue;
                try {
                    done.add(Integer.parseInt(parts[2].trim()) + "|"
                            + Double.parseDouble(parts[1].trim()));
                } catch (NumberFormatException ignore) {}
            }
        }
        return done;
    }

    /** Run one (size, noise) cell: returns {condLines, adjLine}. */
    private static String[] runCell(int size, double noise, int count) throws InitException {
        final double dev = DEV_SCALE * noise;
        final Random rng = ThreadLocalRandom.current();
        final StringBuilder cond = new StringBuilder();

        try { appendCondRow(cond, size, "Hilbert", noise, hilbertMatrix(size, dev, rng)); }
        catch (Exception ignore) {}

        final Stat[][] stats = new Stat[3][3];
        final int[][]  losses = new int [3][3];
        for (int c = 0; c < 3; ++c)
            for (int m = 0; m < 3; ++m) stats[c][m] = new Stat();
        final int[] adjHisto = new int[HISTO_BINS];

        int actualCount = 0;
        for (int r = 0; r < count; ++r) {
            try {
                final int[][] sInt = createIntMatrix(size, rng, ELEMENT_RANGE);
                final Matrix adjPrecise = asMatrix(sInt).adjugate();
                final Matrix mNoisy = addNoise(sInt, dev, rng);
                appendCondRow(cond, size, "Random", noise, mNoisy);
                accumulate(sInt, dev, mNoisy, adjPrecise, stats, losses, adjHisto);
                ++actualCount;
            } catch (Exception ignore) {}
        }

        final StringBuilder adj = new StringBuilder();
        appendAdjRow(adj, noise, size, actualCount, stats, losses, adjHisto);
        return new String[]{cond.toString(), adj.toString()};
    }

    // 31-bin histogram for the normalized Adj error, range [-3, +3] step 0.2.
    private static final int    HISTO_DIVIDS = 5;
    private static final double HISTO_DEVS   = 3.0;
    private static final int    HISTO_HALF   = (int) Math.round(HISTO_DIVIDS * HISTO_DEVS);
    private static final int    HISTO_BINS   = 2 * HISTO_HALF + 1;

    private static void accumHisto(double v, int[] histo) {
        if (v < -HISTO_DEVS || v > HISTO_DEVS) return;
        histo[HISTO_HALF + (int) Math.round(v * HISTO_DIVIDS)]++;
    }

    private static final String CONDITION_HEADER =
        "Size\tType\tNoise\tCondition Number\tDeterminant Value"
      + "\tDeterminant Uncertainty\tDeterminant Precision\tRun Time\n";

    private static final String ADJ_HEADER = buildAdjHeader();
    private static String buildAdjHeader() {
        final String[] CAT = {"Adj", "Fwd", "Rnd"};
        final String[] MET = {"Unc", "Val", "Norm"};
        final String[] FLD = {"Deviation", "Mean", "Minimum", "Maximum", "Count", "Loss"};
        final StringBuilder sb = new StringBuilder("Type\tNoise\tSize\tCount");
        for (String c : CAT) for (String m : MET) for (String f : FLD)
            sb.append('\t').append(c).append(' ').append(m).append(' ').append(f);
        // Adj-norm histogram bin centers: -3.0, -2.8, …, 2.8, 3.0
        for (int i = -HISTO_HALF; i <= HISTO_HALF; ++i)
            sb.append('\t').append((double) i / HISTO_DIVIDS);
        return sb.append('\n').toString();
    }

    private static void appendCondRow(StringBuilder f, int size, String type,
                                       double noise, Matrix m) throws InitException {
        final long t0 = System.nanoTime();
        double cond;
        try { cond = computeCond(m); } catch (Exception e) { cond = Double.NaN; }
        final VarDbl det = m.determ();
        final double val = det.value(), unc = det.uncertainty();
        final double prec = val != 0 ? unc / Math.abs(val) : Double.POSITIVE_INFINITY;
        final double rt = (System.nanoTime() - t0) / 1e9;
        f.append(size).append('\t').append(type).append('\t').append(noise)
         .append('\t').append(cond).append('\t').append(val)
         .append('\t').append(unc).append('\t').append(prec)
         .append('\t').append(rt).append('\n');
    }

    private static void appendAdjRow(StringBuilder f, double noise, int size,
                                      int count, Stat[][] stats, int[][] losses,
                                      int[] adjHisto) {
        f.append("Gaussian\t").append(noise).append('\t').append(size).append('\t').append(count);
        for (int c = 0; c < 3; ++c) {
            for (int m = 0; m < 3; ++m) {
                final Stat s = stats[c][m];
                f.append('\t').append(s.dev()).append('\t').append(s.avg())
                 .append('\t').append(s.min()).append('\t').append(s.max())
                 .append('\t').append(s.count()).append('\t').append(losses[c][m]);
            }
        }
        // 31 bin probabilities for the normalized adj error.
        int total = 0;
        for (int b : adjHisto) total += b;
        for (int i = 0; i < HISTO_BINS; ++i) {
            f.append('\t');
            if (total > 0) f.append((double) adjHisto[i] / total);
        }
        f.append('\n');
    }

    private static double computeCond(Matrix m) {
        final int n = m.size;
        final double[][] data = new double[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) data[i][j] = m.sValue[i * n + j].value();
        final RealMatrix rm = MatrixUtils.createRealMatrix(data);
        return new SingularValueDecomposition(rm).getConditionNumber();
    }

    /** Feed one noisy matrix's Adj/Fwd/Rnd cell errors into the stat grid.
     *  The Adj normalization uses an UNCONDITIONAL predictor — Formula (6.9)
     *  applied to a (sInt means, dev uncertainty) matrix — so the recorded
     *  `Adj Norm` reflects "how well does Formula (6.9) predict the actual
     *  noise-induced spread of adj(M)?". Fwd and Rnd still use the
     *  sample-conditional uncertainty propagation through M·adj. */
    private static void accumulate(int[][] sInt, double dev,
                                    Matrix mNoisy, Matrix adjPrecise,
                                    Stat[][] stats, int[][] losses,
                                    int[] adjHisto) throws InitException {
        final int n = mNoisy.size;

        // Unconditional Adj predictor: Formula (6.9) at sInt means with dev unc.
        final Matrix preciseWithDev = new Matrix(n);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                preciseWithDev.set(r, c, new VarDbl((double) sInt[r][c], dev));
        final Matrix adjPredicted = preciseWithDev.adjugate();

        final Matrix adjNoisy = mNoisy.adjugate();
        final VarDbl detNoisy = mNoisy.determ();
        final double detVal = detNoisy.value(), detVar = detNoisy.variance();
        final Matrix ssId  = mNoisy.multiply(adjNoisy);
        final Matrix ssIdL = adjNoisy.multiply(mNoisy);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                final int idx = i * n + j;
                final VarDbl adjCell = adjNoisy.sValue[idx];
                final VarDbl preCell = adjPrecise.sValue[idx];
                final VarDbl idCell  = ssId.sValue[idx];
                final VarDbl idLCell = ssIdL.sValue[idx];

                // Adj: unconditional predictor → ratio ≈ 1 if Formula (6.9) is faithful.
                final double adjPredUnc = adjPredicted.sValue[idx].uncertainty();
                final double adjErr     = adjCell.value() - preCell.value();
                feedStats(stats[0], losses[0], adjPredUnc, adjErr, adjPredUnc);
                if (adjPredUnc > 0) accumHisto(adjErr / adjPredUnc, adjHisto);

                // Fwd: (M·adj − det·I)[i,j] — rounding error of an identity; the
                // sample-conditional propagated unc is the right denominator here.
                final double idVal = idCell.value(), idUnc = idCell.uncertainty();
                final double fwdUnc = (i == j) ? Math.sqrt(idUnc * idUnc + detVar) : idUnc;
                feedStats(stats[1], losses[1], fwdUnc,
                          (i == j) ? idVal - detVal : idVal, fwdUnc);

                // Rnd: (M·adj − adj·M)[i,j] — same reasoning as Fwd.
                final double idLUnc = idLCell.uncertainty();
                final double rndUnc = Math.sqrt(idUnc * idUnc + idLUnc * idLUnc);
                feedStats(stats[2], losses[2], rndUnc, idVal - idLCell.value(), rndUnc);
            }
        }
    }

    private static void feedStats(Stat[] row, int[] lossesRow,
                                   double unc, double val, double normUnc) {
        row[0].accum(unc, 0);
        row[1].accum(val, 0);
        if (normUnc > 0) row[2].accum(val / normUnc, 0);
        else             lossesRow[2]++;
    }
}
