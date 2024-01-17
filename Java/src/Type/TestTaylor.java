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

    void dumpHeader(final FileWriter fw, final String func) throws IOException {
        fw.write("NoiseType\tNoise\tX\t");
        fw.write(func);
        fw.write("\tDeviation\tUncertainty\tMean\tBias\tLeak");
        for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i) {
            fw.write(String.format("\t%.1f", ((double) i) / DIVIDS));  
        }
        fw.write("\n");
    }

    void dumpResult(final FileWriter fw, double x, double y, 
                    boolean gaussian, double noise, double leak,
                    final Stat stat, final Histogram histo) throws IOException {
        fw.write(String.format("%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", 
                    gaussian? "Gaussian" : "Uniform", noise, x, y, 
                    stat.dev(), var.uncertainty(), 
                    stat.avg() - y, var.value() - y, leak / SAMPLES));
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
    }

    @Test
    public void testPowerExpansion() {
        final double[] sExp = new double[] {-2, -1.5, -1, -0.5, 0.5, 1.5};
        final double[] sDev = new double[] {0.2, 0.195};
        final double[][] ssVar = new double[sExp.length][];
        try (
            final FileWriter fw = new FileWriter("./Java/Output/PowerExpansion.txt")) {
            fw.write("2n\tMomentum\t");
            for (int j = 0; j < sExp.length; ++j) {
                ssVar[j] = new double[sDev.length];
                for (int i = 0; i < sDev.length; ++i) {
                    fw.write(String.format("Exp=%.2f, Dev=%.2e\tVariance\tValue\t", 
                            sExp[j], sDev[i]));
                    ssVar[j][i] = 1;
                }
            }
            fw.write("\n");
            for (int n = 2; n < Momentum.maxN*2; n += 2) {
                fw.write(String.format("%d\t%e\t", n, Momentum.factor(n, BINDING)));
                for (int j = 0; j < sExp.length; ++j) {
                    final double[] sTaylor = Taylor.power(sExp[j]);
                    for (int i = 0; i < sDev.length; ++i) {
                        ssVar[j][i] *= sDev[i] * sDev[i];
                        final double value = Momentum.factor(n, BINDING) * ssVar[j][i];
                        double variance = 0;
                        for (int k = 1; k < n; ++k) {
                            variance += sTaylor[k] * sTaylor[n-k] * 
                                Momentum.factor(n, BINDING) * ssVar[j][i];
                        }
                        for (int k = 2; k < n; k += 2) {
                            variance -= sTaylor[k] * sTaylor[n-k] * 
                                Momentum.factor(k, BINDING) * Momentum.factor(n - k, BINDING) * ssVar[j][i];
                        }
                        fw.write(String.format("%e\t%e\t%e\t", ssVar[j][i], variance, value));
                    }
                }
                fw.write("\n");
            }
            fw.close();
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
            assertEquals(1+1E-2, var.value(), 2E-9);
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
        final double[] sNoise = new double[] {0.2, 0.195, 0.1, 0.05, 0.02, 0.01};
        try (final FileWriter fw = new FileWriter("./Java/Output/PowVar.txt", !gaussian)) {
           if (gaussian) 
                dumpHeader(fw, "Pow");
            for (double dev: sNoise) {
                for (double p: sPower) {
                    try {
                        final double[] sTaylor = Taylor.power(p);
                        sTaylor[0] = 1;
                        var = new VarDbl(1, dev);
                        var.taylor(String.format("pow^(%f)", p), sTaylor, true, true, BINDING);

                        Stat stat = new Stat();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.pow(1 + d*dev, p);
                            assertEquals(1, Math.pow(1 + d*dev, -p)*res, 4E-16);
                            stat.accum(res);
                        }

                        double leak = 0;
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        final double rdev = stat.dev();
                        final double ravg = stat.avg();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.pow(1 + d*dev, p);
                            if (!histo.accum((res - ravg)/rdev)) {
                                leak += 1;
                            }
                        }

                        dumpResult(fw, p, 1, 
                                    gaussian, dev, leak, stat, histo);
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
            }
        } catch (IOException e) {
            fail(e.getMessage());
        } 
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
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(0), 2E-9);
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
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(-1), 2E-9);
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
        final double[] sNoise = new double[] {0.2, 0.1, 0.05, 0.02, 0.01};
        try (final FileWriter fw = new FileWriter("./Java/Output/ExpVar.txt", !gaussian)) {
            if (gaussian)
                dumpHeader(fw, "Sin");

            for (double noise: sNoise) {
                for (double exp: sExp) {
                    try {
                        final double[] sTaylor = Taylor.exp();
                        sTaylor[0] = Math.exp(exp);
                        var = new VarDbl(exp, noise);
                        var.taylor(String.format("exp", exp, noise), sTaylor, false, true, BINDING);

                        Stat stat = new Stat();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.exp(exp + d*noise);
                            stat.accum(res);
                        }

                        double leak = 0;
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        final double rdev = stat.dev();
                        final double ravg = stat.avg();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.exp(exp + d*noise);
                            if (!histo.accum((res - ravg)/rdev)) {
                                leak += 1;
                            }
                        }
                        dumpResult(fw, exp, Math.exp(exp),  
                                    gaussian, noise, leak, stat, histo);
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
            }
        } catch (IOException e) {
            fail(e.getMessage());
        } 
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
            assertEquals(Math.log(1) -1E-2*1/2 -1E-4*3/4 -1E-6*15/6 -1E-8*105/8, var.value(), 2E-9);
            assertEquals(1E-2 +1E-4*5/2 +1E-6*32/3 + 1E-8*65, var.variance() /Math.exp(0)/Math.exp(0), 7E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.log(2);
        init(2, 0.1);
        try {
            var.taylor("log", sTaylor, true, false, 5);
            assertEquals(Math.log(2) -1E-2/4*1/2 -1E-4/16*3/4 -1E-6/64*15/6 -1E-8/256*105/8, var.value(), 2E-9);
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
        final double[] sBase = new double[] {
            0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5, 0.75, 1, 2, 5, 10};
        final double[] sNoise = new double[] {0.2, 0.1, 0.05, 0.02, 0.01};
        try (final FileWriter fw = new FileWriter("./Java/Output/LogVar.txt", !gaussian)) {
            if (gaussian)
                dumpHeader(fw, "Log");

            for (double noise: sNoise) {
                for (double base: sBase) {
                    final double[] sTaylor = Taylor.log();
                    sTaylor[0] = Math.log(base);
                    try {
                        var = new VarDbl(base, noise);
                        var.taylor("log", sTaylor, true, false, BINDING);

                        Stat stat = new Stat();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double x = base + d*noise;
                            if (x <= 0)
                                throw new IllegalArgumentException();
                            final double res = Math.log(x);
                            stat.accum(res);
                        }

                        double leak = 0;
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        final double rdev = stat.dev();
                        final double ravg = stat.avg();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.log(base + d*noise);
                            if (!histo.accum((res - ravg)/rdev)) {
                                leak += 1;
                            }
                        }
                        dumpResult(fw, base, Math.log(base), 
                                    gaussian, noise, leak, stat, histo);
                    } catch (ValueException | UncertaintyException | IllegalArgumentException e) {
                        continue;
                    }
                }
            }
        } catch (IOException e) {
            fail(e.getMessage());
        } 
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
        final double[] sRad = new double[] {
            -1.0, -1.0/12*11, -1.0/6*5, -1.0/4*3, -1.0/3*2, -1.0/12*7, 
            -1.0/2, -1.0/12*5, -1.0/3, -1.0/4, -1.0/6, -1.0/12,
            0, 1.0/12, 1.0/6, 1.0/4, 1.0/3, 1.0/12*5,
            1.0/2, 1.0/12*7, 1.0/3*2, 1.0/4*3, 1.0/6*5, 1.0/12*11, 1.0};
        final double[] sNoise = new double[] {0.2, 0.1, 0.05, 0.02, 0.01};
        try (final FileWriter fw = new FileWriter("./Java/Output/SinVar.txt", !gaussian)) {
            if (gaussian) 
                dumpHeader(fw, "Sin");

            for (double dev: sNoise) {
                for (double rad: sRad) {
                    try {
                        final double[] sTaylor = Taylor.sin(rad*Math.PI);
                        var = new VarDbl(rad*Math.PI, dev);
                        var.taylor("sin", sTaylor, false, false, BINDING);

                        Stat stat = new Stat();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.sin(rad*Math.PI + d*dev);
                            stat.accum(res);
                        }

                        double leak = 0;
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        final double rdev = stat.dev();
                        final double ravg = stat.avg();
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? sGass[i] : sUnif[i];
                            final double res = Math.sin(rad*Math.PI + d*dev);
                            if (!histo.accum((res - ravg)/rdev)) {
                                leak += 1;
                            }
                        }
                        dumpResult(fw, rad, Math.sin(rad*Math.PI), 
                                    gaussian, dev, leak, stat, histo);
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
            }
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }
}


