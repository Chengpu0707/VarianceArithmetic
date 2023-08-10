package Type;

import static org.junit.Assert.assertEquals;
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
    public void testPower() {
        double[] sTaylor;

        sTaylor = Taylor.power(0);
        assertEquals(0, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("(1~0.1)^0", sTaylor, true, true, 5);
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
            var.taylor("(1~0.1)^1", sTaylor, true, true, 5);
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
            var.taylor("(10~0.1)^1", sTaylor, true, true, 5);
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
            var.taylor("(1~0.1)^2", sTaylor, true, true, 5);
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
            var.taylor("(1~0.1)^(1/2)", sTaylor, true, true, 5);
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
            var.taylor("(1~0.1)^(-1)", sTaylor, true, true, 5);
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
            var.taylor("(1~0.1)^(-2)", sTaylor, true, true, 5);
            assertEquals(1 +1E-2*3 +1E-4*15 +1E-6*105, var.value(), 2E-4);
            assertEquals(1E-2*4 +1E-4*66 + 1E-6*960, var.variance(), 2E-3);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void TestPowerDump() {
        System.out.println(System.getProperty("user.dir"));
        final int SAMPLES = 10000, BINDING = 5, DIVIDS = 2;
        final double[] sPower = new double[] {
            -2, -1.75, -5.0/3, -1.5, -4.0/3, -1.25, 
            -1, -0.75, -2.0/3, -0.5, -1.0/3, -0.25, -0.1, -0.01 -1E-3, -1E-6, 
            1E-6, 1E-3, 0.01, 0.1, 0.25, 1.0/3, 0.5, 2.0/3, 0.75, 
            1, 1.25, 4.0/3, 1.5, 5.0/3, 1.75, 
            2, 2.25, 7.0/3, 2.5, 8.0/3, 2.75, 3};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
        Random rand = new Random();
        final StringBuilder sb = new StringBuilder();
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/Power.txt")) {
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
                        var.taylor(String.format("(1~%f)^%f", dev, p), sTaylor, true, true, BINDING);
                        final double unc = var.uncertainty();
                        final double mode = (1 + Math.sqrt(1 - 4*(p-1)*dev*dev))/2;

                        double leak = 0;
                        Stat stat = new Stat();
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = rand.nextGaussian();
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
            var.taylor("e^(0~0.1)", sTaylor, false, true, 5);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(0), 2E-9);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(0)/Math.exp(0), 7E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.exp(1);
        init(1, 0.1);
        try {
            var.taylor("e^(1~0.1)", sTaylor, false, true, 5);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(1), 1E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(1)/Math.exp(1), 6E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor[0] = Math.exp(-1);
        init(-1, 0.1);
        try {
            var.taylor("e^(1~0.1)", sTaylor, false, true, 5);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(-1), 2E-9);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(-1)/Math.exp(-1), 6E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void TestExpDump() {
        System.out.println(System.getProperty("user.dir"));
        final int SAMPLES = 10000, BINDING = 5, DIVIDS = 2;
        final double[] sExp = new double[] {
            -100, -50, -20, -10, -5, -2, -1, -0.5, -0.2, -0.1,
            0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
        Random rand = new Random();
        final StringBuilder sb = new StringBuilder();
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/Exp.txt")) {
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

            for (double p: sExp) {
                fw.write(String.format("%f\t", p));
                sb.setLength(0);
                for (double dev: sDev) {
                    try {
                        final double[] sTaylor = Taylor.exp();
                        sTaylor[0] = Math.exp(p);
                        final VarDbl var = new VarDbl(p, dev * dev);
                        var.taylor(String.format("e^(%f~%f)^", p, dev), sTaylor, false, true, BINDING);
                        final double unc = var.uncertainty();
                        final double mode = p - dev*dev;

                        double leak = 0;
                        Stat stat = new Stat();
                        Histogram histo = new Histogram(BINDING, DIVIDS);
                        for (int i = 0; i < SAMPLES; ++i) {
                            final double d = rand.nextGaussian();
                            if (BINDING <= Math.abs(d)) {
                                leak += 1;
                                continue;
                            }
                            final double res = Math.exp(p + d*dev);
                            stat.accum(res/Math.exp(p));
                            histo.accum((res - Math.exp(p))/unc);
                        }
                        fw.write(String.format("%g\t%g\t%g\t%g\t%g\t%g\t", 
                                 leak / SAMPLES, stat.avg(), var.value() /Math.exp(p), stat.dev(), 
                                 unc /Math.exp(p) /Math.exp(p), mode));
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


