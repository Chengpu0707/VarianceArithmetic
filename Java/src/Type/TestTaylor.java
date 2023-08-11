package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;

import java.util.Random;

import org.junit.Test;

import Stats.Histogram;
import Stats.Stat;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;


public class TestTaylor {
    final static double TOLERANCE = 3E-16;
    VarDbl var;

    private void init(double value, double dev) {
        try {
            var = new VarDbl(value, dev*dev);
        } catch (ValueException | UncertaintyException e) {
            fail();
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
        final double[] sTaylor = Taylor.power(-2);
        final double[] sDev = new double[] {0.2, 0.199, 0.198, 0.197, 0.196, 0.195};
        final double bounding = 5;
        final double[] sVar = new double[sDev.length];
        try (
            final FileWriter fw = new FileWriter("./PowerExpansion.txt")) {
            fw.write("2n\tMomentum\t");
            for (int i = 0; i < sDev.length; ++i) {
                fw.write(String.format("%e\t%.3f Value\t%.3f Variance\t", sDev[i], sDev[i], sDev[i]));
                sVar[i] = 1;
            }
            fw.write("\n");
            for (int n = 2; n < Momentum.maxN*2; n += 2) {
                fw.write(String.format("%d\t%e\t", n, Momentum.factor(n, bounding)));
                for (int i = 0; i < sDev.length; ++i) {
                    sVar[i] *= sDev[i] * sDev[i];
                    final double value = Momentum.factor(n, bounding) * sVar[i];
                    double variance = 0;
                    for (int j = 1; j < n; ++j) {
                        variance += sTaylor[j] * sTaylor[n-j] * Momentum.factor(n, bounding) * sVar[i];
                    }
                    for (int j = 2; j < n; j += 2) {
                        variance -= sTaylor[j] * sTaylor[n-j] * Momentum.factor(j, bounding) * Momentum.factor(n - j, bounding) * sVar[i];
                    }
                    fw.write(String.format("%e\t%e\t%e\t", sVar[i], value, variance));
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
    public void TestPowerDump() {
        final boolean gaussian = true;
        final int SAMPLES = 10000, BINDING = 5, DIVIDS = 2;
        final double[] sPower = new double[] {
            -2, -1.75, -5.0/3, -1.5, -4.0/3, -1.25, 
            -1, -0.75, -2.0/3, -0.5, -1.0/3, -0.25, -0.1, -0.01 -1E-3, -1E-6, 
            1E-6, 1E-3, 0.01, 0.1, 0.25, 1.0/3, 0.5, 2.0/3, 0.75, 
            1, 1.25, 4.0/3, 1.5, 5.0/3, 1.75, 
            2, 2.25, 7.0/3, 2.5, 8.0/3, 2.75, 3};
        final double[] sDev = new double[] {0.2, 0.195, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002};
        Random rand = new Random();
        final StringBuilder sb = new StringBuilder();
        try (
            final FileWriter fw = new FileWriter(gaussian? "./PowerGN.txt" : "./PowerUN.txt")) {
            fw.write("Dev\t");
            for (double dev: sDev) {
                fw.write(String.format("%g\t\t\t\t\t\t", dev));
                sb.append(String.format("%g\t", dev));
                for (int i = -BINDING*DIVIDS; i < BINDING*DIVIDS; ++i) {
                    sb.append("\t");  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            sb.setLength(0);
            fw.write("exponent\t");
            for (double dev: sDev) {
                fw.write("Leak\tMean\tmean\tDev\tdev\tmode\t");
                for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i) {
                    sb.append(String.format("%g\t", ((double) i) / DIVIDS));  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            for (double p: sPower) {
                fw.write(String.format("%f\t", p));
                sb.setLength(0);
                for (double dev: sDev) {
                    try {
                        final double[] sTaylor = Taylor.power(p);
                        sTaylor[0] = 1;
                        final VarDbl var = new VarDbl(1, dev * dev);
                        var.taylor(String.format("pow^(%f)", p), sTaylor, true, true, BINDING);
                        final double unc = var.uncertainty();
                        final double mode = (1 + Math.sqrt(1 - 4*(p-1)*dev*dev))/2;

                        double leak = 0;
                        Stat stat = new Stat();
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? rand.nextGaussian() : (rand.nextDouble() - 0.5) * Math.sqrt(12);
                            if (BINDING <= Math.abs(d)) {
                                leak += 1;
                                continue;
                            }
                            final double res = Math.pow(1 + d*dev, p);
                            assertEquals(1, Math.pow(1 + d*dev, -p)*res, 4E-16);
                            stat.accum(res);
                            histo.accum((res - var.value())/unc);
                        }
                        fw.write(String.format("%g\t%g\t%g\t%g\t%g\t%g\t", 
                                 leak / SAMPLES, stat.avg(), var.value(), stat.dev(), unc, mode));
                        final double[] sHisto = histo.histo();
                        for (int i = 0; i < sHisto.length; ++i) {
                            sb.append(String.format("%g\t", sHisto[i]));
                        }
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
                fw.write(sb.toString());
                fw.write("\n");
            }
            fw.close();
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
    public void TestExpDump() {
        final boolean gaussian = true;
        final int SAMPLES = 10000, BINDING = 5, DIVIDS = 2;
        final double[] sExp = new double[] {
            -100, -50, -20, -10, -5, -2, -1, -0.5, -0.2, -0.1,
            0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
        Random rand = new Random();
        final StringBuilder sb = new StringBuilder();
        try (
            final FileWriter fw = new FileWriter(gaussian? "./ExpGN.txt" : "./ExpUN.txt")) {
            fw.write("Dev\t");
            for (double dev: sDev) {
                fw.write(String.format("%g\t\t\t\t\t\t", dev));
                sb.append(String.format("%g\t", dev));
                for (int i = -BINDING*DIVIDS; i < BINDING*DIVIDS; ++i) {
                    sb.append("\t");  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            sb.setLength(0);
            fw.write("exponent\t");
            for (double dev: sDev) {
                fw.write("Leak\tMean\tmean\tDev\tdev\tmode\t");
                for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i) {
                    sb.append(String.format("%g\t", ((double) i) / DIVIDS));  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            for (double exp: sExp) {
                fw.write(String.format("%f\t", exp));
                sb.setLength(0);
                for (double dev: sDev) {
                    try {
                        final double[] sTaylor = Taylor.exp();
                        sTaylor[0] = Math.exp(exp);
                        final VarDbl var = new VarDbl(exp, dev * dev);
                        var.taylor(String.format("exp", exp, dev), sTaylor, false, true, BINDING);
                        final double unc = var.uncertainty();
                        final double mode = exp - dev*dev;

                        double leak = 0;
                        Stat stat = new Stat();
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? rand.nextGaussian() : (rand.nextDouble() - 0.5) * Math.sqrt(12);
                            if (BINDING <= Math.abs(d)) {
                                leak += 1;
                                continue;
                            }
                            final double res = Math.exp(exp + d*dev);
                            stat.accum(res);
                            histo.accum((res - Math.exp(exp))/unc);
                        }
                        fw.write(String.format("%g\t%g\t%g\t%g\t%g\t%g\t", 
                                 leak / SAMPLES, stat.avg(), var.value(), stat.dev(), unc, mode));
                        final double[] sHisto = histo.histo();
                        for (int i = 0; i < sHisto.length; ++i) {
                            sb.append(String.format("%g\t", sHisto[i]));
                        }
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
                fw.write(sb.toString());
                fw.write("\n");
            }
            fw.close();
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
    public void TestLogDump() {
        final boolean gaussian = true;
        final int SAMPLES = 10000, BINDING = 5, DIVIDS = 2;
        final double[] sBase = new double[] {
            0.5, 0.75, 0.875, 
            1, 2, 5, 10};
        final double[] sDev = new double[] {0.195, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
        Random rand = new Random();
        final StringBuilder sb = new StringBuilder();
        try (
            final FileWriter fw = new FileWriter(gaussian? "./LogGN.txt" : "./LogUN.txt")) {
            fw.write("Dev\t");
            for (double dev: sDev) {
                fw.write(String.format("%g\t\t\t\t\t\t", dev));
                sb.append(String.format("%g\t", dev));
                for (int i = -BINDING*DIVIDS; i < BINDING*DIVIDS; ++i) {
                    sb.append("\t");  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            sb.setLength(0);
            fw.write("exponent\t");
            for (double dev: sDev) {
                fw.write("Leak\tMean\tmean\tDev\tdev\tmode\t");
                for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i) {
                    sb.append(String.format("%g\t", ((double) i) / DIVIDS));  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            for (double base: sBase) {
                fw.write(String.format("%f\t", base));
                sb.setLength(0);
                for (double dev: sDev) {
                    try {
                        final double[] sTaylor = Taylor.log();
                        sTaylor[0] = Math.log(base);
                        final VarDbl var = new VarDbl(base, dev * dev);
                        var.taylor("log", sTaylor, true, false, BINDING);
                        final double unc = var.uncertainty();
                        final double mode = (base + Math.sqrt(base*base + 4*dev*dev))/2;

                        double leak = 0;
                        Stat stat = new Stat();
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? rand.nextGaussian() : (rand.nextDouble() - 0.5) * Math.sqrt(12);
                            if (BINDING <= Math.abs(d)) {
                                leak += 1;
                                continue;
                            }
                            if (base + d*dev <= 0) {
                                leak += 1;
                                continue;
                            }
                            final double res = Math.log(base + d*dev);
                            stat.accum(res);
                            histo.accum(res - Math.log(base));
                        }
                        fw.write(String.format("%g\t%g\t%g\t%g\t%g\t%g\t", 
                                 leak / SAMPLES, stat.avg(), var.value(), stat.dev(), unc, mode));
                        final double[] sHisto = histo.histo();
                        for (int i = 0; i < sHisto.length; ++i) {
                            sb.append(String.format("%g\t", sHisto[i]));
                        }
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
                fw.write(sb.toString());
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }

    @Test
    public void TestSin() {
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
    public void TestSinDump() {
        final boolean gaussian = true;
        final int SAMPLES = 10000, BINDING = 5, DIVIDS = 2;
        final double[] sRad = new double[] {
            -1.0, -1.0/12*11, -1.0/6*5, -1.0/4*3, -1.0/3*2, -1.0/12*7, 
            -1.0/2, -1.0/12*5, -1.0/3, -1.0/4, -1.0/6, -1.0/12,
            0, 1.0/12, 1.0/6, 1.0/4, 1.0/3, 1.0/12*5,
            1.0/2, 1.0/12*7, 1.0/3*2, 1.0/4*3, 1.0/6*5, 1.0/12*11, 1.0};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
        Random rand = new Random();
        final StringBuilder sb = new StringBuilder();
        try (
            final FileWriter fw = new FileWriter(gaussian? "./SinGN.txt" : "./SinUN.txt")) {
            fw.write("Dev\t");
            for (double dev: sDev) {
                fw.write(String.format("%g\t\t\t\t\t\t", dev));
                sb.append(String.format("%g\t", dev));
                for (int i = -BINDING*DIVIDS; i < BINDING*DIVIDS; ++i) {
                    sb.append("\t");  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            sb.setLength(0);
            fw.write("sin\t");
            for (double dev: sDev) {
                fw.write("Leak\tMean\tmean\tDev\tdev\tMode\t");
                for (int i = -BINDING*DIVIDS; i <= BINDING*DIVIDS; ++i) {
                    sb.append(String.format("%g\t", ((double) i) / DIVIDS /8));  
                }
            }
            fw.write(sb.toString());
            fw.write("\n");

            for (double rad: sRad) {
                fw.write(String.format("%f\t", rad));
                sb.setLength(0);
                for (double dev: sDev) {
                    try {
                        final double[] sTaylor = Taylor.sin(rad*Math.PI);
                        final VarDbl var = new VarDbl(rad*Math.PI, dev * dev);
                        var.taylor("sin", sTaylor, false, false, BINDING);
                        final double unc = var.uncertainty();
                        final double mode = var.value() - Math.sin(rad*Math.PI);

                        double leak = 0;
                        Stat stat = new Stat();
                        Histogram histo = new Histogram(BINDING/8.0, DIVIDS*8);
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = gaussian? rand.nextGaussian() : (rand.nextDouble() - 0.5) * Math.sqrt(12);
                            if (BINDING <= Math.abs(d)) {
                                leak += 1;
                                continue;
                            }
                            final double res = Math.sin(rad*Math.PI + d*dev);
                            stat.accum(res);
                            histo.accum(res - var.value());
                        }
                        fw.write(String.format("%g\t%g\t%g\t%g\t%g\t%g\t", 
                                 leak / SAMPLES, stat.avg(), var.value(), stat.dev(), unc, mode));
                        final double[] sHisto = histo.histo();
                        for (int i = 0; i < sHisto.length; ++i) {
                            sb.append(String.format("%g\t", sHisto[i]));
                        }
                    } catch (ValueException | UncertaintyException e) {
                        fail(e.getMessage());
                    }
                }
                fw.write(sb.toString());
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }
}


