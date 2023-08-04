package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

import org.junit.Test;

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

        sTaylor = Taylor.powerTaylor(0);
        assertEquals(0, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("(1~0.1)^0", sTaylor, true, 5);
            assertEquals(1, var.value(), 0);
            assertEquals(0, var.variance(), 0);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.powerTaylor(1);
        assertEquals(1, sTaylor[1], TOLERANCE);
        assertEquals(0, sTaylor[2], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("(1~0.1)^1", sTaylor, true, 5);
            assertEquals(1, var.value(), 0);
            assertEquals(0.01, var.variance(), 2E-6);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.powerTaylor(2);
        assertEquals(2, sTaylor[1], TOLERANCE);
        assertEquals(1, sTaylor[2], TOLERANCE);
        assertEquals(0, sTaylor[3], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("(1~0.1)^2", sTaylor, true, 5);
            assertEquals(1+1E-2, var.value(), 2E-9);
            assertEquals(1E-2*4+1E-4*2, var.variance(), 5E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.powerTaylor(0.5);
        assertEquals(1.0/2, sTaylor[1], TOLERANCE);
        assertEquals(-1.0/8, sTaylor[2], TOLERANCE);
        assertEquals(1.0/16, sTaylor[3], TOLERANCE);
        assertEquals(-5.0/128, sTaylor[4], TOLERANCE);
        assertEquals(7.0/256, sTaylor[5], TOLERANCE);
        assertEquals(-21.0/1024, sTaylor[6], TOLERANCE);
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("(1~0.1)^(1/2)", sTaylor, true, 5);
            assertEquals(1 -1E-2/8 -1E-4*15/128 -1E-6*315/1024, var.value(), 2E-7);
            assertEquals(1E-2/4 +1E-4*7/32 +1E-6*75/128, var.variance(), 3E-7);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }

        sTaylor = Taylor.powerTaylor(-1);
        for (int i = 1; i < sTaylor.length; ++i) {
            assertEquals(((i % 2) == 1)? -1 : 1, sTaylor[i], TOLERANCE);
        }
        sTaylor[0] = 1;
        init(1, 0.1);
        try {
            var.taylor("(1~0.1)^(-1)", sTaylor, true, 5);
            assertEquals(1 +1E-2 +1E-4*3 +1E-6*15, var.value(), 2E-6);
            assertEquals(1E-2 +1E-4*8 + 1E-6*69, var.variance(), 1E-4);
        } catch (ValueException | UncertaintyException e) {
            fail();
        }
    }

    @Test
    public void TestPowerDump() {
        System.out.println(System.getProperty("user.dir"));
        final int SAMPLES = 10000, BINDING = 5;
        final double[] sPower = new double[] {
            -2, -1.75, -5.0/3, -1.5, -4.0/3, -1.25, 
            -1, -0.75, -2.0/3, -0.5, -1.0/3, -0.25, -0.1, -0.01 -1E-3, -1E-6, 
            1E-6, 1E-3, 0.01, 0.1, 0.25, 1.0/3, 0.5, 2.0/3, 0.75, 
            1, 1.25, 4.0/3, 1.5, 5.0/3, 1.75, 
            2, 2.25, 7.0/3, 2.5, 8.0/3, 2.75, 3};
        final double[] sDev = new double[] {0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
        Random rand = new Random();
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/Power.txt")) {
            fw.write("Dev\t");
            for (double dev: sDev) {
                fw.write(String.format("%f\t\t\t\t\t\t", dev));
            }
            fw.write("\n");
            fw.write("exponent\t");
            for (double dev: sDev) {
                fw.write("Leak\tMean\tmean\tDev\tdev\tDev/dev\t");
            }
            fw.write("\n");

            for (double p: sPower) {
                fw.write(String.format("%f\t", p));
                for (double dev: sDev) {
                    double leak = 0;
                    Stat stat = new Stat();
                    for (int i = 0; i < SAMPLES; ++i) {
                        final double d = rand.nextGaussian();
                        if (BINDING <= Math.abs(d)) {
                            leak += 1;
                            continue;
                        }
                        stat.accum(Math.pow(1 + d*dev, p));
                    }
                    try {
                        final double[] sTaylor = Taylor.powerTaylor(p);
                        sTaylor[0] = 1;
                        final VarDbl var = new VarDbl(1, dev * dev);
                        var.taylor(String.format("(1~%f)^%f", dev, p), sTaylor, true, BINDING);
                        final double unc = var.uncertainty();
                        fw.write(String.format("%g\t%g\t%g\t%g\t%g\t%g\t", 
                                 leak / SAMPLES, stat.avg(), var.value(), stat.dev(), unc, (unc == 0)? 1 : stat.dev()/unc));
                    } catch (ValueException | UncertaintyException e) {
                        fail();
                    }
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }

}


