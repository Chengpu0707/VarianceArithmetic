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
            f.write("X\tEdge\tBias\tValue\tUncertainty\tException\n");
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
                                res = Taylor.pow(new VarDbl(1, dx), x);
                                val = 1;
                                break;
                            case "Exp":
                                res = Taylor.exp(new VarDbl(x, dx));
                                val = Math.exp(x);
                                break;
                            case "Log":
                                res = Taylor.log(new VarDbl(x, dx));
                                val = Math.log(x);
                                break;
                            case "Sin":
                                res = Taylor.sin(new VarDbl(x, dx));
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
                f.write(String.format("%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%s\n", 
                    x, edge, res.value() - val, res.value(), res.uncertainty(), except));}
        } catch (IOException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testPow() throws InitException {
        try {
            Taylor.pow(new VarDbl(1, 0.19929), -1.75, "./Java/Output/pow_1_0.19929_-1.75.txt");
        } catch (NotStableException ex) {
        } catch (Throwable e) {
            fail(String.format("(1+-0.19929)^-1.75 throws %s", e));
        }

        final double[] sX = new double[701];
        for (int i = -60; i <= 80; i += 1) {
            sX[i + 60] = i/20.;
        }
        searchEdge("Pow", sX, new int[]{18000, 22000, 100000}, true);
    }

    @Test
    public void testSin() throws InitException {
        final int DIVIDS = 64;
        final double[] sX = new double[DIVIDS * 2 + 1];
        for (int i = -DIVIDS; i <= DIVIDS; i += 1) {
            sX[i + DIVIDS] = Math.PI * i /DIVIDS;
        }
        searchEdge("Sin", sX, new int[]{9000, 20000, 10000}, false);
    }

    @Test
    public void testLog() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        final double[] sX = new double[]{
            1./64, 1./32, 1./16, 1./8, 1./4, 1./2, 1, 2, 4, 8, 16, 32, 64
        };
        for (double x: sX) {
            try {
                final VarDbl res = Taylor.log(new VarDbl(x, x * 0.20086));
                assertEquals(res.uncertainty(), 0.2130506, 1e-6);
            } catch (Throwable ex) {
                fail(String.format("Log(%f) fails for %s", x, ex.getMessage()));
            }
            try {
                Taylor.log(new VarDbl(x, x * 0.20087));
                fail(String.format("Log(%f) not fails", x));
            } catch (NotMonotonicException ex) {
            } catch (IOException e) {
                fail(String.format("Log(%f) fails: %s", x, e));
            }
        }
    }

    @Test
    public void testExp() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        final double[] sX = new double[]{
            -50, -20, -10, -5, -2, 0, 2, 5, 10, 20, 50
        };
        for (double x: sX) {
            try {
                final VarDbl res = Taylor.exp(new VarDbl(x, 19.864));
                assertEquals(res.uncertainty() / res.value(), 1681.767, 1e-3);
            } catch (Throwable ex) {
                fail(String.format("Exp(%f) fails for %s", x, ex.getMessage()));
            }
            try {
                Taylor.exp(new VarDbl(x, 19.865));
                fail(String.format("Exp(%f) not fails", x));
            } catch (NotMonotonicException ex) {
            } catch (IOException e) {
                fail(String.format("Exp(%f) fails for %s", x, e.getMessage()));
            }
        }
    }
}
