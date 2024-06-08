package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.Random;

import Stats.Histogram;
import Stats.Stat;


public class TestTaylor {
    static final int HIST_RANGE = 3;
    static final int HIST_DIVIDS = 5;
    static final int SAMPLES = 10000;

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
        } catch (InitException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testEnvVar() {
        assertNull(System.getProperty("workspaceFolder"));
        assertNull(System.getenv("workspaceFolder"));
        System.out.println(System.getProperty("user.dir"));
    }

    void dumpTest(final String test, boolean gaussian, final double[] sX, final double[] sDev) {
        try (final FileWriter fw = new FileWriter(String.format("./Java/Output/%sVar.txt", test), !gaussian)) {
            if (gaussian) {
                fw.write("NoiseType\tNoise\tX\t");
                fw.write(test);
                fw.write("\tError Deviation\tError Minimum\tError Maximum\tValue Deviation\tUncertainty\tMean\tBias\tLeak");
                for (int i = -HIST_RANGE * HIST_DIVIDS; 
                        i <=HIST_RANGE * HIST_DIVIDS; ++i)
                    fw.write(String.format("\t%.1f", ((double) i) / HIST_DIVIDS));  
                fw.write("\n");
            }
            for (double dev: sDev) {
                for (double x: sX) {
                    double expect;
                    try {
                        final VarDbl[] sTaylor;
                        if (test.equals("pow")) {
                            sTaylor = Taylor.power(x);
                            sTaylor[0] = new VarDbl(1);
                            var = new VarDbl(1, dev);
                            var = var.power(x);
                        } else if (test == "exp") {
                            sTaylor = Taylor.exp();
                            sTaylor[0] = new VarDbl(Math.exp(x));  
                            var = new VarDbl(x, dev);
                            var = var.taylor(String.format("exp(%s)", var), sTaylor, false, true);
                        } else if (test == "log") {
                            sTaylor = Taylor.log();
                            sTaylor[0] = new VarDbl(Math.log(x));
                            var = new VarDbl(x, dev);
                            var = var.taylor(String.format("log(%s)", var), sTaylor, true, false);
                        } else if (test == "sin") {
                            sTaylor = Taylor.sin(x*Math.PI);
                            sTaylor[0] = new VarDbl(Math.sin(x*Math.PI));
                            var = new VarDbl(x*Math.PI, dev);
                            var = var.taylor(String.format("sin(%s)", var), sTaylor, false, false);
                        } else {
                            throw new IllegalArgumentException(String.format("Unkonwn test %s", test));
                        }
                        expect = sTaylor[0].value(); 
                    } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException e) {
                        System.out.println(e.getMessage());
                        continue;
                    }

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
                    Histogram histo = new Histogram(HIST_RANGE, HIST_DIVIDS);
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
                        stat.accum((res - expect)/var.uncertainty());
                        if (!histo.accum((res - vavg)/vdev)) {
                            leak += 1;
                        }
                    }

                    fw.write(String.format("%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", 
                                gaussian? "Gaussian" : "Uniform", dev, x, expect, 
                                stat.dev(), stat.min(), stat.max(), vdev, var.uncertainty(), 
                                vavg - expect, var.value() - expect, leak / SAMPLES));
                    final double[] sHisto = histo.histo();
                    if (sHisto == null) {
                        for (int i = -HIST_RANGE * HIST_DIVIDS; 
                                i <= HIST_RANGE * HIST_DIVIDS ; ++i) {
                            fw.write("\t");
                        }
                    } else {
                        for (int i = 0; i < sHisto.length; ++i) {
                            fw.write(String.format("\t%g", sHisto[i]));
                        }
                    }
                    fw.write("\n");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            fail(e.getMessage());
        }  
    }

    @Test
    public void testPower_0() {
        try {
            final VarDbl[] sTaylor = Taylor.power(0);
            assertEquals(0, sTaylor[1].value(), 0);
            assertNull(sTaylor[2]);

            init(1, 0.1);
            var = var.power(0);
            assertEquals(1, var.value(), 0);
            assertEquals(0, var.uncertainty(), 0);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException | IOException e) {
            fail();
        }
    }

    @Test
    public void testPower_1() {
        try {
            final VarDbl[] sTaylor = Taylor.power(1);
            assertEquals(1, sTaylor[1].value(), 0);
            assertEquals(0, sTaylor[2].value(), 0);
            assertNull(sTaylor[3]);

            init(1, 0.1);
            var = var.power(1);
            assertEquals(1, var.value(), 0);
            assertEquals(0.1, var.uncertainty(), 2E-6);

            init(10, 0.1);
            var = var.power(1);
            assertEquals(10, var.value(), 0);
            assertEquals(0.01, var.variance(), 2E-6);

        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException | IOException e) {
            fail();
        }
    }

    @Test
    public void testPower_2() {
        try {
            final VarDbl[] sTaylor = Taylor.power(2);
            assertEquals(2, sTaylor[1].value(), 0);
            assertEquals(1, sTaylor[2].value(), 0);
            assertEquals(0, sTaylor[3].value(), 0);
            assertNull(sTaylor[4]);

            init(1, 0.1);
            var = var.power(2, "./Java/Output/Power2.txt");
            assertEquals(1 +1E-2, var.value(), 2E-7);
            assertEquals(1E-2*4 +1E-4*2, var.variance(), 1E-6);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException | IOException e) {
            fail();
        }
    }

    @Test
    public void testPower_half() {
        try {
            final VarDbl[] sTaylor = Taylor.power(0.5);
            assertEquals(1.0/2, sTaylor[1].value(), 0);
            assertEquals(-1.0/8, sTaylor[2].value(), 0);
            assertEquals(1.0/16, sTaylor[3].value(), 0);
            assertEquals(-5.0/128, sTaylor[4].value(), 0);
            assertEquals(7.0/256, sTaylor[5].value(), 0);
            assertEquals(-21.0/1024, sTaylor[6].value(), 0);
            
            init(1, 0.1);
            var = var.power(0.5);
            assertEquals(1 -1E-2/8 -1E-4*15/128 -1E-6*315/1024, var.value(), 2E-7);
            assertEquals(1E-2/4 +1E-4*7/32 +1E-6*75/128, var.variance(), 3E-7);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException | IOException e) {
            fail();
        }
    }

    @Test
    public void testPower_inv() {
        try {
            final VarDbl[] sTaylor = Taylor.power(-1);
            for (int i = 1; i < sTaylor.length; ++i) {
                assertEquals(((i % 2) == 1)? -1 : 1, sTaylor[i].value(), 0);
            }
            
            init(1, 0.1);
            var = var.power(-1);
            assertEquals(1 +1E-2 +1E-4*3 +1E-6*15, var.value(), 2E-6);
            assertEquals(1E-2 +1E-4*8 + 1E-6*69, var.variance(), 1E-4);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException | IOException e) {
            fail();
        }
    }

    @Test
    public void testPower_inv2() {
        try {
            final VarDbl[] sTaylor = Taylor.power(-2);
            for (int i = 1; i < sTaylor.length; ++i) {
                assertEquals(((i % 2) == 1)? -(i+1) : (i+1), sTaylor[i].value(), sTaylor[i].uncertainty());
            }
            
            init(1, 0.1);
            var = var.power(-2);
            assertEquals(1 +1E-2*3 +1E-4*15 +1E-6*105, var.value(), 2E-4);
            assertEquals(1E-2*4 +1E-4*66 + 1E-6*960, var.variance(), 2E-3);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException | IOException e) {
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
        VarDbl[] sTaylor;

        try {
            sTaylor = Taylor.exp();
            assertEquals(0, sTaylor[0].value(), 0);
            assertEquals(1.0/1, sTaylor[1].value(), 0);
            assertEquals(1.0/2, sTaylor[2].value(), 0);
            assertEquals(1.0/6, sTaylor[3].value(), 0);
            assertEquals(1.0/24, sTaylor[4].value(), 0);
            assertEquals(1.0/120, sTaylor[5].value(), 0);

            sTaylor[0] = new VarDbl(Math.exp(0));
            init(0, 0.1);
            var = var.taylor("exp", sTaylor, false, true);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(0), 2E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(0)/Math.exp(0), 7E-7);

            sTaylor[0] = new VarDbl(Math.exp(1));
            init(1, 0.1);
            var = var.taylor("exp", sTaylor, false, true);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(1), 1E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6, var.variance() /Math.exp(1)/Math.exp(1), 6E-7);

            sTaylor[0] = new VarDbl(Math.exp(-1));
            init(-1, 0.1);
            var = var.taylor("exp", sTaylor, false, true);
            assertEquals(1 +1E-2/2 +1E-4/8 +1E-6/48, var.value() /Math.exp(-1), 2E-7);
            assertEquals(1E-2 +1E-4*3/2 + 1E-6*7/6+1E-8*5/8, var.variance() /Math.exp(-1)/Math.exp(-1), 6E-7);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException e) {
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
        final double[] sDev = new double[] {0.2, 0.1, 0.01, 0.001, 1E-4, 1E-5, 1E-6};
        dumpTest("exp", gaussian, sExp, sDev);
    }

    @Test
    public void testLog() {
        try {
            final VarDbl[] sTaylor = Taylor.log();
            assertEquals(0, sTaylor[0].value(), 0);
            assertEquals(1.0/1, sTaylor[1].value(), 0);
            assertEquals(-1.0/2, sTaylor[2].value(), 0);
            assertEquals(1.0/3, sTaylor[3].value(), 0);
            assertEquals(-1.0/4, sTaylor[4].value(), 0);
            assertEquals(1.0/5, sTaylor[5].value(), 0);

            sTaylor[0] = new VarDbl(Math.log(1));
            init(1, 0.1);
            var = var.taylor("log", sTaylor, true, false);
            assertEquals(Math.log(1) -1E-2*1/2 -1E-4*3/4 -1E-6*15/6 -1E-8*105/8, var.value(), 2E-7);
            assertEquals(1E-2 +1E-4*5/2 +1E-6*32/3 + 1E-8*65, var.variance() /Math.exp(0)/Math.exp(0), 7E-7);

            sTaylor[0] = new VarDbl(Math.log(2));
            init(2, 0.1);
            var = var.taylor("log", sTaylor, true, false);
            assertEquals(Math.log(2) -1E-2/4*1/2 -1E-4/16*3/4 -1E-6/64*15/6 -1E-8/256*105/8, var.value(), 2E-7);
            assertEquals(1E-2/4 +1E-4/16*5/2 +1E-6/64*32/3 + 1E-8/256*65, var.variance() /Math.exp(0)/Math.exp(0), 5E-7);

            sTaylor[0] = new VarDbl(Math.log(0.5));
            init(0.5, 0.1);
            var = var.taylor("log", sTaylor, true, false);
            assertEquals(Math.log(0.5) -4E-2*1/2 -16E-4*3/4 -64E-6*15/6 -256E-8*105/8, var.value(), 2E-4);
            assertEquals(4E-2 +16E-4*5/2 +64E-6*32/3 + 256E-8*65, var.variance(), 1E-4);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException e) {
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
        final double[] sDev = new double[] {0.2, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6};
        dumpTest("log", gaussian, sX, sDev);
    }

    @Test
    public void testSin() {
        VarDbl[] sTaylor;
        try {
            sTaylor = Taylor.sin(Math.PI / 6);
            assertEquals(0.5, sTaylor[0].value(), Math.ulp(sTaylor[0].value()));
            assertEquals(Math.sqrt(3)/2, sTaylor[1].value(), Math.ulp(sTaylor[1].value()));
            assertEquals(-0.5 /2, sTaylor[2].value(), Math.ulp(sTaylor[2].value()));
            assertEquals(-Math.sqrt(3)/2 /6, sTaylor[3].value(), Math.ulp(sTaylor[3].value()));
            assertEquals(0.5 /24, sTaylor[4].value(), Math.ulp(sTaylor[4].value()));
            assertEquals(Math.sqrt(3)/2 /120, sTaylor[5].value(), Math.ulp(sTaylor[5].value()));
            init(Math.PI / 6, 0.1);
            var = var.taylor("sin", sTaylor, false, false);
            assertEquals(0.5 -1E-2*0.5/2 +1E-4*0.5/24 -1E-6*0.5/720 +1E-8*0.5/264320, var.value(), 5E-5);
            assertEquals(1E-2*3/4 -1E-4*5/8 +1E-6*23/96, var.variance(), 2E-6);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException e) {
            fail();
        }

        try {
            sTaylor = Taylor.sin(0);
            assertEquals(0, sTaylor[0].value(), 0);
            assertEquals(1, sTaylor[1].value(), 0);
            assertEquals(0, sTaylor[2].value(), 0);
            assertEquals(-1.0 /6, sTaylor[3].value(), 0);
            assertEquals(0, sTaylor[4].value(), 0);
            assertEquals(1.0 /120, sTaylor[5].value(), 0);
            init(0, 0.1);
            var = var.taylor("sin", sTaylor, false, false);
            assertEquals(0, var.value(), 0);
            assertEquals(1E-2 -1E-4 +1E-6*13/24, var.variance(), 2E-6);
        } catch (InitException | DivergentException | NotReliableException | NotMonotonicException | NotStableException e) {
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
        final double[] sDev = new double[] {0.2, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6};
        dumpTest("sin", gaussian, sX, sDev);
    }


    @Test
    public void testExpansion() {
        try {
            new VarDbl(1, 0.2).power(-1, "./Java/Output/Power_1_0.2_-1.txt");
            fail();
        } catch (InitException | DivergentException | NotReliableException | NotStableException | IOException e) {
            fail(e.getMessage());
        } catch (NotMonotonicException e) {
            assertEquals(102, e.order);
            try (BufferedReader br = new BufferedReader(new FileReader("./Java/Output/Power_1_0.2_-1.txt"))) {
                String line, last = null;
                int cnt = 0;
                while ((line = br.readLine()) != null) {
                    ++cnt;
                    last = line;
                }
                assertEquals(58, cnt);
                assertEquals("NotMonotonicException", last);
            } catch (IOException ex) {
                fail(ex.getMessage());
            }
        }
        
        try {
            VarDbl res = new VarDbl(2, 0.2).power(-1, "./Java/Output/Power_2_0.2_-1.txt");
            try (BufferedReader br = new BufferedReader(new FileReader("./Java/Output/Power_2_0.2_-1.txt"))) {
                String line, last = null;
                int cnt = 0;
                while ((line = br.readLine()) != null) {
                    ++cnt;
                    last = line;
                }
                assertEquals(26, cnt);
                final String[] sLast = last.split("\t");
                assertEquals(res.value(), Double.valueOf(sLast[0]), 1e-7);          
                assertEquals(res.variance(), Double.valueOf(sLast[1]) + Double.valueOf(sLast[2]), 1e-8);
            } catch (IOException e) {
                fail(e.getMessage());
            }
        } catch (InitException | DivergentException | NotMonotonicException | NotReliableException | NotStableException | IOException e) {
            fail(e.getMessage());
        }
    }
}


