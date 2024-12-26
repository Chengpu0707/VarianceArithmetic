package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.FileWriter;
import java.io.IOException;


public class TestConvergence {
    private void searchEdge(final String type, final double[] sX, final int[] sSearch, boolean ignoreNoBound) 
            throws InitException {
        try (FileWriter f = new FileWriter(String.format("./Java/Output/%sEdge.txt", type))) {
            f.write("Edge Value\tEdge Uncertainty\tBias\tUncertainty\tException\n");
            assertNotEquals(0, sX.length);
            assertEquals(3, sSearch.length);
            assertTrue(sSearch[0] < sSearch[1]);
            assertTrue(0 < sSearch[2]);
            for (double x: sX) {
                int i = sSearch[0], j = sSearch[1];
                double val = 0, edge = 0;
                VarDbl res = null;
                String except = null;
                while (i + 1 < j) {
                    final int k = (i + j)/2;
                    final double dx = ((double) k)/sSearch[2];
                    try {
                        switch (type) {
                            case "Pow":
                                res = new VarDbl(1, dx).pow(x);
                                val = 1;
                                break;
                            case "Exp":
                                res = new VarDbl(x, dx).exp();
                                val = Math.exp(x);
                                break;
                            case "Log":
                                res = new VarDbl(x, dx).log();
                                val = Math.log(x);
                                break;
                            case "Sin":
                                res = new VarDbl(x, dx).sin();
                                val = Math.sin(x);
                                break;
                            default:
                                throw new IllegalArgumentException(String.format("Unknown type %s", type));
                        }
                        edge = dx;
                        i = k;
                    } catch (Taylor1dException ex) {
                        except = ex.getClass().getName();
                        j = k;
                    }
                }
                assertNotEquals(i, sSearch[0]);
                if (ignoreNoBound && (j == sSearch[1]))
                    continue; 
                assertNotEquals(j, sSearch[1]);
                f.write(String.format("%e\t%e\t%e\t%e\t%s\n", 
                    x, edge, res.value() - val, res.uncertainty(), except));}
        } catch (IOException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testPow() throws InitException {
        final double[] sX = new double[701];
        for (int i = -300; i <= 400; i += 1) {
            sX[i + 300] = i/100.;
        }
        searchEdge("Pow", sX, new int[]{18000, 22000, 100000}, true);
    }

    @Test
    public void testSin() throws InitException {
        final int DIVIDS = 32;
        final double[] sX = new double[DIVIDS * 2 + 1];
        for (int i = -DIVIDS; i <= DIVIDS; i += 1) {
            sX[i + DIVIDS] = Math.PI * i /DIVIDS;
        }
        searchEdge("Sin", sX, new int[]{9000, 15000, 10000}, false);
    }

    @Test
    public void testLog() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        final double[] sX = new double[]{
            1./64, 1./32, 1./16, 1./8, 1./4, 1./2, 1, 2, 4, 8, 16, 32, 64
        };
        for (double x: sX) {
            try {
                final VarDbl res = new VarDbl(x, x * 0.200887).log();
                assertEquals(res.uncertainty(), 0.213083, 1e-6);
            } catch (Throwable ex) {
                fail(String.format("Log(%f) fails for %s", x, ex.getMessage()));
            }
            try {
                new VarDbl(x, x * 0.200888).log();
                fail(String.format("Log(%f) not fails", x));
            } catch (NotMonotonicException ex) {

            }
        }
    }

    @Test
    public void testExp() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        final double[] sX = new double[]{
            -100, -50, -20, -10, -5, -2, 0, 2, 5, 10, 20, 50, 100
        };
        for (double x: sX) {
            try {
                final VarDbl res = new VarDbl(x, 19.864).exp();
                assertEquals(res.uncertainty() / res.value(), 1680.377, 1e-3);
            } catch (Throwable ex) {
                fail(String.format("Log(%f) fails for %s", x, ex.getMessage()));
            }
            try {
                new VarDbl(x, 19.865).exp();
                fail(String.format("Log(%f) not fails", x));
            } catch (NotMonotonicException ex) {
            }
        }
    }
}
