package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.FileWriter;
import java.io.IOException;

import java.util.Random;

import Stats.Histogram;
import Stats.Stat;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;


public class TestTaylor {
    static final int SAMPLES = 10000;
    static final int BINDING = 5, DIVIDS = 5;
    static final double TOLERANCE = 3E-16;

    static final double[] sGass = new double[SAMPLES];
    static final double[] sUnif = new double[SAMPLES];
    static {
        final double delta = 1.0/(SAMPLES - 1);
        for (var i = 0; i < SAMPLES; ++i) {
            Random rand = new Random();
            sGass[i] = rand.nextGaussian();
            sUnif[i] = (delta * i - 0.5) * Math.sqrt(12);
        }
    }
    VarDbl var;

    private void init(double value, double dev) {
        try {
            var = new VarDbl(value, dev);
        } catch (ValueException | UncertaintyException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testEnvVar() {
        assertNull(System.getProperty("workspaceFolder"));
        assertNull(System.getenv("workspaceFolder"));
        System.out.println(System.getProperty("user.dir"));
    }

    @Test
    public void testPowerExpansion() {
        final double[] sExp = new double[] {-2, -1.5, -1, -0.5, 0.5, 1.5};
        final double[] sDev = new double[] {0.2, 0.195};
        try (
            final FileWriter fw = new FileWriter("./Java/Output/PowerExpansion.txt")) {
            fw.write("2n\tMomentum\tExponent\tInput Uncertainty\tTerm\tVariance\n");
            fw.write("\n");
             for (int n = 2; n < Momentum.maxN*2; n += 2) {
                final double factor = Momentum.factor(n, BINDING);
                for (int j = 0; j < sExp.length; ++j) {
                    final double[] sTaylor = Taylor.power(sExp[j]);
                    for (int i = 0; i < sDev.length; ++i) {
                        double variance = 0;
                        for (int k = 1; k < n; ++k) {
                            variance += sTaylor[k] * sTaylor[n-k] * factor;
                        }
                        for (int k = 2; k < n; k += 2) {
                            variance -= sTaylor[k] * sTaylor[n-k] * 
                                Momentum.factor(k, BINDING) * Momentum.factor(n - k, BINDING);
                        }
                        fw.write(String.format("%d\t%e\t%e\t%e\t%e\t%e\n", 
                                 n, factor, sExp[j], sDev[i], variance, variance * Math.pow(sDev[i], n)));
                    }
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }


    void dumpTest(final String test, boolean gaussian, final double[] sX, final double[] sDev) {
        try (final FileWriter fw = new FileWriter(String.format("./Java/Output/%sVar.txt", test), !gaussian)) {
            if (gaussian) {
                fw.write("NoiseType\tNoise\tX\t");
                fw.write(test);
                fw.write("\tError Deviation\tError Minimum\tError Maximum\tValue Deviation\tUncertainty\tMean\tBias\tLeak");
                for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i)
                    fw.write(String.format("\t%.1f", ((double) i) / DIVIDS));  
                fw.write("\n");
            }
            for (double dev: sDev) {
                for (double x: sX) {
                    try {
                        final double[] sTaylor;
                        if (test.equals("pow")) {
                            sTaylor = Taylor.power(x);
                            sTaylor[0] = 1;
                            var = new VarDbl(1, dev);
                            var.taylor(test, sTaylor, true, true, BINDING);
                        } else if (test == "exp") {
                            sTaylor = Taylor.exp();
                            sTaylor[0] = Math.exp(x);  
                            var = new VarDbl(x, dev);
                            var.taylor(test, sTaylor, false, true, BINDING);
                        } else if (test == "log") {
                            sTaylor = Taylor.log();
                            sTaylor[0] = Math.log(x);
                            var = new VarDbl(x, dev);
                            var.taylor(test, sTaylor, true, false, BINDING);
                        } else if (test == "sin") {
                            sTaylor = Taylor.sin(x*Math.PI);
                            sTaylor[0] = Math.sin(x*Math.PI);
                            var = new VarDbl(x*Math.PI, dev);
                            var.taylor(test, sTaylor, false, false, BINDING);
                        } else
                            throw new IllegalArgumentException(String.format("Unkonwn test %s", test));
 
                        Stat stat = new Stat();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res;
                            if (test.equals("pow"))
                                res = Math.pow(1 + d*dev, x);
                            else if (test.equals("exp"))
                                res = Math.exp(x + d*dev);
                            else if (test.equals("log"))
                                res = Math.log(x + d*dev);
                            else if (test == "sin")
                                res = Math.sin(x*Math.PI + d*dev);
                            else
                                throw new IllegalArgumentException(String.format("Unkonwn test %s", test));
                            stat.accum(res);
                        }
 
                         double leak = 0;
                         Histogram histo = new Histogram(BINDING, DIVIDS);
                         final double vdev = stat.dev();
                         final double vavg = stat.avg();
                         stat.clear();
                         for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res;
                            if (test.equals("pow"))
                                res = Math.pow(1 + d*dev, x);
                            else if (test.equals("exp"))
                                res = Math.exp(x + d*dev);
                            else if (test.equals("log"))
                                res = Math.log(x + d*dev);
                            else if (test == "sin")
                                res = Math.sin(x*Math.PI + d*dev);
                            else
                                throw new IllegalArgumentException(String.format("Unkonwn test %s", test));
                            stat.accum((res - sTaylor[0])/var.uncertainty());
                            if (!histo.accum((res - vavg)/vdev)) {
                                leak += 1;
                            }
                        }
 
                        fw.write(String.format("%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", 
                                    gaussian? "Gaussian" : "Uniform", dev, x, sTaylor[0], 
                                    stat.dev(), stat.min(), stat.max(), vdev, var.uncertainty(), 
                                    vavg - sTaylor[0], var.value() - sTaylor[0], leak / SAMPLES));
                        final double[] sHisto = histo.histo();
                        if (sHisto == null) {
                            for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i) {
                                fw.write("\t");
                            }
                        } else {
                            for (int i = 0; i < sHisto.length; ++i) {
                                fw.write(String.format("\t%g", sHisto[i]));
                            }
                        }
                        fw.write("\n");
                    } catch (ValueException | UncertaintyException e) {
                        if (!test.equals("log"))
                            fail(e.getMessage());
                    }
                 }
             }
         } catch (IOException e) {
             fail(e.getMessage());
         }  
    }

    @Test
    public void testPower() {
        double[] sTaylor;

        sTaylor = Taylor.power(0);
        assertEquals(0, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("pow^0", sTaylor, true, true, 5);
            assertEquals(1, var.value(), 0);
            assertEquals(0, var.variance(), 0);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.power(1);
        assertEquals(1, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("pow^1", sTaylor, true, true, 5);
            assertEquals(1, var.value(), 0);
            assertEquals(0.01, var.variance(), 2E-6);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.power(1);
        assertEquals(1, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        sTaylor[0] = 10;
        init(10, 0.1);
        try {
            var.taylor("pow^1", sTaylor, true, true, 5);
            assertEquals(10, var.value(), 0);
            assertEquals(0.01, var.variance(), 2E-6);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.power(2);
        assertEquals(2, sTaylor[1], TOLERANCE);
        assertEquals(1, sTaylor[2], TOLERANCE);
        assertEquals(0, sTaylor[3], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("pow^2", sTaylor, true, true, 5);
            assertEquals(1+1E-2, var.value(), 2E-7);
            assertEquals(1E-2*4+1E-4*2, var.variance(), 5E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.power(0.5);
        assertEquals(1.0/2, sTaylor[1], TOLERANCE);
        assertEquals(-1.0/8, sTaylor[2], TOLERANCE);
        assertEquals(1.0/16, sTaylor[3], TOLERANCE);
        assertEquals(-5.0/128, sTaylor[4], TOLERANCE);
        assertEquals(7.0/256, sTaylor[5], TOLERANCE);
        assertEquals(-21.0/1024, sTaylor[6], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("pow^(1/2)", sTaylor, true, true, 5);
            assertEquals(1 -1E-2/8 -1E-4*15/128 -1E-6*315/1024, var.value(), 2E-7);
            assertEquals(1E-2/4 +1E-4*7/32 +1E-6*75/128, var.variance(), 3E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.power(-1);
        for (int i = 1; i < sTaylor.length; ++i) {
            assertEquals(((i % 2) == 1)? -1 : 1, sTaylor[i], TOLERANCE);
        }
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("pow^(-1)", sTaylor, true, true, 5);
            assertEquals(1 +1E-2 +1E-4*3 +1E-6*15, var.value(), 2E-6);
            assertEquals(1E-2 +1E-4*8 + 1E-6*69, var.variance(), 1E-4);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.power(-2);
        for (int i = 1; i < sTaylor.length; ++i) {
            assertEquals(((i % 2) == 1)? -(i+1) : (i+1), sTaylor[i], TOLERANCE);
        }
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("pow^(-2)", sTaylor, true, true, 5);
            assertEquals(1 +1E-2*3 +1E-4*15 +1E-6*105, var.value(), 2E-4);
            assertEquals(1E-2*4 +1E-4*66 + 1E-6*960, var.variance(), 2E-3);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testPowDump() {
        testPowDump(true);
        testPowDump(false);
    }
    private void testPowDump(boolean gaussian) {
        final double[] sPower = new double[] {
            -2, -1.75, -5.0/3, -1.5, -4.0/3, -1.25, -1.001,
            -1, -0.999, -0.75, -2.0/3, -0.5, -1.0/3, -0.25, -0.1, -0.01 -1E-3, -1E-6, 
            1E-6, 1E-3, 0.01, 0.1, 0.25, 1.0/3, 0.5, 2.0/3, 0.75, 
            1, 1.25, 4.0/3, 1.5, 5.0/3, 1.75,
            2, 2.25, 7.0/3, 2.5, 8.0/3, 2.75, 3};
        final double[] sDev = new double[] {0.2, 0.195, 0.1, 0.05, 0.02, 0.01};
        dumpTest("pow", gaussian, sPower, sDev);
    }

    @Test
    public void testExp() {
        final double[] sTaylor = Taylor.exp();
        assertEquals(0, sTaylor[0], TOLERANCE);
        assertEquals(1.0/1, sTaylor[1], TOLERANCE);
        assertEquals(1.0/2, sTaylor[2], TOLERANCE);
        assertEquals(1.0/6, sTaylor[3], TOLERANCE);
        assertEquals(1.0/24, sTaylor[4], TOLERANCE);
        assertEquals(1.0/120, sTaylor[5], TOLERANCE);

        sTaylor[0] = Math.exp(0);
        init(0, 0.1);
        try {
            var.taylor("exp", sTaylor, false, true, 5);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(0), 2E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(0)/Math.exp(0), 7E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.exp(1);
        init(1, 0.1);
        try {
            var.taylor("exp", sTaylor, false, true, 5);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(1), 1E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(1)/Math.exp(1), 6E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.exp(-1);
        init(-1, 0.1);
        try {
            var.taylor("exp", sTaylor, false, true, 5);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(-1), 2E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6+1E-8*5/8, var.variance() /Math.exp(-1)/Math.exp(-1), 6E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testExpDump() {
        testExpDump(true);
        testExpDump(false);
    }
    private void testExpDump(boolean gaussian) {
        final double[] sExp = new double[] {-10, -5, -2, -1, 0, 1, 2, 5, 10};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01};
        dumpTest("exp", gaussian, sExp, sDev);
    }

    @Test
    public void testLog() {
        final double[] sTaylor = Taylor.log();
        assertEquals(0, sTaylor[0], TOLERANCE);
        assertEquals(1.0/1, sTaylor[1], TOLERANCE);
        assertEquals(-1.0/2, sTaylor[2], TOLERANCE);
        assertEquals(1.0/3, sTaylor[3], TOLERANCE);
        assertEquals(-1.0/4, sTaylor[4], TOLERANCE);
        assertEquals(1.0/5, sTaylor[5], TOLERANCE);

        sTaylor[0] = Math.log(1);
        init(1, 0.1);
        try {
            var.taylor("log", sTaylor, true, false, 5);
            assertEquals(Math.log(1) -1E-2*1/2 -1E-4*3/4 -1E-6*15/6 -1E-8*105/8, var.value(), 2E-7);
            assertEquals(1E-2 +1E-4*5/2 +1E-6*32/3 + 1E-8*65, var.variance() /Math.exp(0)/Math.exp(0), 7E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.log(2);
        init(2, 0.1);
        try {
            var.taylor("log", sTaylor, true, false, 5);
            assertEquals(Math.log(2) -1E-2/4*1/2 -1E-4/16*3/4 -1E-6/64*15/6 -1E-8/256*105/8, var.value(), 2E-7);
            assertEquals(1E-2/4 +1E-4/16*5/2 +1E-6/64*32/3 + 1E-8/256*65, var.variance() /Math.exp(0)/Math.exp(0), 2E-8);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.log(0.5);
        init(0.5, 0.1);
        try {
            var.taylor("log", sTaylor, true, false, 5);
            assertEquals(Math.log(0.5) -4E-2*1/2 -16E-4*3/4 -64E-6*15/6 -256E-8*105/8, var.value(), 2E-4);
            assertEquals(4E-2 +16E-4*5/2 +64E-6*32/3 + 256E-8*65, var.variance(), 1E-4);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testLogDump() {
        testLogDump(true);
        testLogDump(false);
    }
    private void testLogDump(boolean gaussian) {
        final double[] sX = new double[] {
            0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5, 0.75, 1, 2, 5, 10};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01};
        dumpTest("log", gaussian, sX, sDev);
    }

    @Test
    public void testSin() {
        double[] sTaylor = Taylor.sin(Math.PI / 6);
        assertEquals(0.5, sTaylor[0], TOLERANCE);
        assertEquals(Math.sqrt(3)/2, sTaylor[1], TOLERANCE);
        assertEquals(-0.5 /2, sTaylor[2], TOLERANCE);
        assertEquals(-Math.sqrt(3)/2 /6, sTaylor[3], TOLERANCE);
        assertEquals(0.5 /24, sTaylor[4], TOLERANCE);
        assertEquals(Math.sqrt(3)/2 /120, sTaylor[5], TOLERANCE);
        init(Math.PI / 6, 0.1);
        try {
            var.taylor("sin", sTaylor, false, false, 5);
            assertEquals(0.5 -1E-2*0.5/2 +1E-4*0.5/24 -1E-6*0.5/720 +1E-8*0.5/264320, var.value(), 5E-5);
            assertEquals(1E-2*3/4 -1E-4*5/8 +1E-6*23/96, var.variance(), 2E-6);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.sin(0);
        assertEquals(0, sTaylor[0], TOLERANCE);
        assertEquals(1, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        assertEquals(-1.0 /6, sTaylor[3], TOLERANCE);
        assertEquals(0, sTaylor[4], TOLERANCE);
        assertEquals(1.0 /120, sTaylor[5], TOLERANCE);
        init(0, 0.1);
        try {
            var.taylor("sin", sTaylor, false, false, 5);
            assertEquals(0, var.value(), 0);
            assertEquals(1E-2 -1E-4 +1E-6*13/24, var.variance(), 2E-6);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void testSinDump() {
        testSinDump(true);
        testSinDump(false);
    }
    private void testSinDump(boolean gaussian) {
        final double[] sX = new double[] {
            -1.0, -1.0/12*11, -1.0/6*5, -1.0/4*3, -1.0/3*2, -1.0/12*7, 
            -1.0/2, -1.0/12*5, -1.0/3, -1.0/4, -1.0/6, -1.0/12,
            0, 1.0/12, 1.0/6, 1.0/4, 1.0/3, 1.0/12*5,
            1.0/2, 1.0/12*7, 1.0/3*2, 1.0/4*3, 1.0/6*5, 1.0/12*11, 1.0};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01};
        dumpTest("sin", gaussian, sX, sDev);
    }
}


