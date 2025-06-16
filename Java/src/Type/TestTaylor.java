package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import org.junit.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

import Stats.Histogram;
import Stats.Stat;



public class TestTaylor {
    static final int HIST_RANGE = 3;
    static final int HIST_DIVIDS = 5;
    static final int SAMPLES = 10000;
    static final double REPR_DELTA = 1e-6;

    VarDbl var;

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
                        if (test.equals("pow")) {
                            var = new VarDbl(1, dev).pow(x);
                            expect = 1;
                        } else if (test == "exp") {
                            var = new VarDbl(x, dev).exp();
                            expect = Math.exp(x);
                        } else if (test == "log") {
                            var = new VarDbl(x, dev).log();
                            expect = Math.log(x);
                        } else if (test == "sin") {
                            var = new VarDbl(x*Math.PI, dev).sin();
                            expect = Math.sin(x);
                        } else {
                            throw new IllegalArgumentException(String.format("Unkonwn test %s", test));
                        }
                    } catch (InitException | NotFiniteException | NotReliableException | NotMonotonicException | NotStableException e) {
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
        } catch (Throwable e) {
            e.printStackTrace();
            fail(e.getMessage());
        }  
    }

    private void validatePow(final double exp, final double x, final double dx, final double valPrec, double varPrec, final String exception) 
            throws IOException { 
        final String dumpFile = String.format("./Java/Output/pow_%s_%s_%s.txt", x, dx, exp);
        try {
            var = new VarDbl(x, dx).pow(exp, dumpFile);
        } catch (InitException e) {
            assertEquals(exception, "InitException");
            return;
        } catch (IllegalArgumentException e) {
            assertEquals(exception, "IllegalArgumentException");
            return;
        } catch (Taylor1dException e) {
            assertEquals(e.getClass().getName(), "Type." + exception);
        }
        final double r = Math.pow(x, exp);
        File f = new File(dumpFile);
        if (!f.exists())
            return;
        VarDbl.DumpResult res = VarDbl.readDumpFile(dumpFile, var, exception);
        final double[] sTaylor = res.s1dTaylor.sDbl;
        assertEquals(r, sTaylor[0], REPR_DELTA);
        if ((exp > 0) && (Math.floor(exp) == Math.ceil(exp))) {
            final int n = (int) exp;
            assertTrue(sTaylor.length < VarDbl.momentum.maxOrder);
            assertEquals(sTaylor[1], n * Math.pow(x, n - 1), REPR_DELTA);
            assertEquals(sTaylor[n], 1, REPR_DELTA);
            for (int i = n + 1; i < sTaylor.length; ++i)
            assertEquals(sTaylor[i], 0, 0);
        } else {
            assertEquals(sTaylor.length, VarDbl.momentum.maxOrder);
            assertEquals(exp/1, sTaylor[1], Math.ulp(sTaylor[1]));
            assertEquals(exp*(exp - 1)/2, sTaylor[2], Math.ulp(sTaylor[2]));
        }
        if (exception != null)
            return;
        if (x == 0) {
            assertTrue(0 <= exp);
            assertTrue(Math.floor(exp) == Math.ceil(exp));
            int n = (int) exp;
            double value = Math.pow(dx, n) * VarDbl.momentum.get(2);
            double variance = 0;
            n *= 2;
            for (int i = 1; i < n; ++i) {
                variance += sTaylor[i] * sTaylor[n - i] * (VarDbl.momentum.get(n) - VarDbl.momentum.get(i) * VarDbl.momentum.get(n - i));
            }
            variance *= Math.pow(dx, n);
            if ((value == 0))
                assertEquals(var.value(), value, valPrec);
            else
                assertEquals(var.value() /value, 1, valPrec);
            if ((variance == 0))
                assertEquals(var.variance(), variance, varPrec);
            else
                assertEquals(var.variance() /variance, 1, varPrec);
        } else {
            final double pr = dx / x;
            final double value = 1 + 
                    Math.pow(pr, 2) * exp*(exp-1)/2. + 
                    Math.pow(pr, 4) * exp*(exp-1)*(exp-2)*(exp-3)/24. +
                    Math.pow(pr, 6) * exp*(exp-1)*(exp-2)*(exp-3)*(exp-4)*(exp-5)/720.;
            final double variance = Math.pow(pr, 2) * exp*exp + 
                    Math.pow(pr, 4) * exp*exp*(exp-1)*(exp-5./3)*3/2 +
                    Math.pow(pr, 6) * exp*exp*(exp-1)*(exp-2)*(exp-2)*(exp-16./7)*7./6;
            if ((r == 0) || (value == 0))
                assertEquals(var.value(), r, valPrec);
            else
                assertEquals(var.value()/r /value, 1, valPrec);
            if ((r == 0) || (variance == 0))
                assertEquals(var.variance(), variance, varPrec);
            else
                assertEquals(var.variance()/r/r /variance, 1, varPrec);
        }
    }

    @Test
    public void testPow_0() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, IOException {
        validatePow(0, 1, 0.1, 0, 0, null);
        validatePow(0, 0.1, 1, 0, 0, null);
        validatePow(0, 10, 1, 0, 0, null);
        validatePow(0, 1, 10, 0, 0, null);
        validatePow(0, 0, 0.1, 0, 0, null);
        validatePow(0, 0, 1, 0, 0, null);
        validatePow(0, 0, 10, 0, 0, null);
        validatePow(0, -1, 0.1, 0, 0, null);
        validatePow(0, -0.1, 1, 0, 0, null);
        validatePow(0, -10, 1, 0, 0, null);
        validatePow(0, -1, 10, 0, 0, null);
    }

    @Test
    public void testPow_1() 
            throws IOException {
        validatePow(1, 1, 0.1, 0, 0, null);
        validatePow(1, 0.1, 1, 0, 0, null);
        validatePow(1, 10, 1, 0, 3e-16, null);
        validatePow(1, 1, 10, 0, 0, null);
        validatePow(1, 0, 0.1, 0, 0, null);
        validatePow(1, 0, 1, 0, 0, null);
        validatePow(1, 0, 10, 0, 0, null);
        validatePow(1, -1, 0.1, 0, 0, null);
        validatePow(1, -0.1, 1, 0, 0, null);
        validatePow(1, -10, 1, 0, 3e-16, null);
        validatePow(1, -1, 10, 0, 0, null);
    }

    @Test
    public void testPow_2() 
            throws IOException {
        validatePow(2, 1, 0.1, 2e-7, 2e-5, null);
        validatePow(2, 0.1, 1, 2e-5, 2e-4, null);
        validatePow(2, 10, 1, 2e-7, 2e-5, null);
        validatePow(2, 1, 10, 2e-5, 2e-4, null);
        validatePow(2, 0, 0.1, 1e-2, 3e-16, null);
        validatePow(2, 0, 1, 2e-5, 0, null);
        validatePow(2, 0, 10, 2e-5, 0, null);
        validatePow(2, -1, 0.1, 2e-7, 2e-5, null);
        validatePow(2, -0.1, 1, 2e-5, 2e-4, null);
        validatePow(2, -10, 1, 2e-7, 2e-5, null);
        validatePow(2, -1, 10, 2e-5, 2e-4, null);
    }

    @Test
    public void testPow_3() 
            throws IOException {
        validatePow(3.01, 1, 0.2023, 7e-6, 5e-5, null);
        validatePow(3.02, 1, 0.2023, 2e-5, 5e-5, null);
        validatePow(3.02, 1, 0.1900, 2e-5, 5e-5, null);
    }

    @Test
    public void testPow_half() 
            throws IOException {
        validatePow(0.5, 1, 0.1, 9e-6, 6e-6, null);
        validatePow(0.5, 0.1, 1, 0, 0, "NotFiniteException");
        validatePow(0.5, 10, 1, 9e-6, 6e-6, null);
        validatePow(0.5, 1, 10, 0, 0, "NotFiniteException");
        validatePow(0.5, 0, 0.1, 0, 0, "NotMonotonicException");
        validatePow(0.5, -1, 0.1, 0, 0, "IllegalArgumentException");
    }

    @Test
    public void testPow_neg1() 
            throws IOException {
        validatePow(-1, 1, 0.1, 3e-4, 7e-4, null);
        validatePow(-1, 10, 1, 3e-4, 7e-4, null);
        validatePow(-1, 0.1, 1, 0, 0, "NotFiniteException");
        validatePow(-1, 1, 0.20003, 0, 0, "NotMonotonicException");
        validatePow(-1, 1, 0.20002, 5e-3, 2e-1, null);
    }

    @Test
    public void testPow_neg2() 
            throws IOException {
        validatePow(-2, 1, 0.1, 2e-3, 4e-3, null);
        validatePow(-2, 10, 1, 2e-3, 4e-3, null);
        validatePow(-2, 0.1, 1, 0, 0, "NotFiniteException");
        validatePow(-2, 1, 0.19906, 0, 0, "NotMonotonicException");
        validatePow(-2, 1, 0.19905, 3e-2, 30, null);
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

     
    private void validateExp(final double x, final double dx, final double precVal, final double precVar, final String exception) 
            throws IOException {
        final String dumpFile = String.format("./Java/Output/exp_%s_%s.txt", x, dx);
        try {
            var = new VarDbl(x, dx).exp(dumpFile);
        } catch (InitException e) {
            assertEquals(exception, "InitException");
        } catch (Taylor1dException e) {
            assertEquals(e.getClass().getName(), "Type." + exception);
        }
        assertEquals(1 + Math.pow(dx, 2)/2 + Math.pow(dx, 4)/8 + Math.pow(dx, 6)/48 + Math.pow(dx, 8)/384, 
                    var.value() /Math.exp(x), precVal);
        assertEquals(Math.pow(dx, 2) + Math.pow(dx, 4) *3/2 + Math.pow(dx, 6) *7/6 + Math.pow(dx, 8) *5/8, 
                    var.variance() /Math.exp(x)/Math.exp(x), precVar);
        VarDbl.DumpResult res = VarDbl.readDumpFile(dumpFile, var, exception);
        final double[] sTaylor = res.s1dTaylor.sDbl;
        assertEquals(VarDbl.momentum.maxOrder, sTaylor.length);
        assertEquals(Math.exp(x), sTaylor[0], REPR_DELTA);
        assertEquals(1.0/1, sTaylor[1], REPR_DELTA);
        assertEquals(1.0/2, sTaylor[2], REPR_DELTA);
        assertEquals(1.0/6, sTaylor[3], REPR_DELTA);
        assertEquals(1.0/24, sTaylor[4], REPR_DELTA);
        assertEquals(1.0/120, sTaylor[5], REPR_DELTA);
    }

    @Test
    public void testExp() 
            throws IOException {
        validateExp(0, 0.1, 2e-7, 7e-7, null);
        validateExp(1, 0.1, 1e-7, 6e-7, null);
        validateExp(-1, 0.1, 2e-7, 6e-7, null);
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


    private void validateLog(final double x, final double dx, final double deltaVal, final double deltaVar, final String exception) 
            throws IOException {
        final String dumpFile = String.format("./Java/Output/log_%s_%s.txt", x, dx);
        try {
            var = new VarDbl(x, dx).log(dumpFile);
        } catch (InitException e) {
            assertEquals(exception, "InitException");
        } catch (Taylor1dException e) {
            assertEquals(e.getClass().getName(), "Type." + exception);
        }
        final double pr = dx / x;
        assertEquals( - Math.pow(pr, 2)/2 - Math.pow(pr, 4)*3/4 - Math.pow(pr, 6)*15/6 - Math.pow(dx, 8)*105/8, 
                    var.value() - Math.log(x), deltaVal);
        assertEquals(Math.pow(pr, 2) + Math.pow(pr, 4) *9/8 + Math.pow(pr, 6) *119/24 + Math.pow(pr, 8) *991/32, 
                    var.variance(), deltaVar);
        VarDbl.DumpResult res = VarDbl.readDumpFile(dumpFile, var, exception);
        final double[] sTaylor = res.s1dTaylor.sDbl;
        assertEquals(VarDbl.momentum.maxOrder, sTaylor.length);
        assertEquals(sTaylor[0], Math.log(x), REPR_DELTA);
        assertEquals(1.0/1, sTaylor[1], REPR_DELTA);
        assertEquals(-1.0/2, sTaylor[2], REPR_DELTA);
        assertEquals(1.0/3, sTaylor[3], REPR_DELTA);
        assertEquals(-1.0/4, sTaylor[4], REPR_DELTA);
        assertEquals(1.0/5, sTaylor[5], REPR_DELTA);
}

    @Test
    public void testLog() 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, InitException, IOException {
        validateLog(1, 0.1, 8e-8, 2e-4, null);
        validateLog(2, 0.1, 2e-7, 9e-6, null);
        validateLog(0.5, 0.1, 5e-5, 3e-3, null);
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


    void validateSin(final double x, final double dx, final double deltaVal, final double deltaVar, final String exception) 
            throws IOException {
        final String dumpFile = String.format("./Java/Output/sin_%s_%s.txt", x, dx);
        try {
            var = new VarDbl(x, dx).sin(dumpFile);
        } catch (InitException e) {
            assertEquals(exception, "InitException");
        } catch (Taylor1dException e) {
            assertEquals(e.getClass().getName(), "Type." + exception);
        }
        final double sin = Math.sin(x);
        final double cos = Math.cos(x);
        assertEquals(sin * (1 - Math.pow(dx, 2)/2 + Math.pow(dx, 4)/8 - Math.pow(dx, 6)/48 + Math.pow(dx, 8)/384), 
                    var.value(), deltaVal);
        final double cos2 = cos*cos;
        assertEquals(Math.pow(dx, 2) *cos2 - Math.pow(dx, 4) *(cos2*3/2 - 1./2) 
                        + Math.pow(dx, 6) *(cos2*7/6 - 1./2) - Math.pow(dx, 8) *(cos2*5/8 - 7./24), 
                    var.variance(), deltaVar);
        VarDbl.DumpResult res = VarDbl.readDumpFile(dumpFile, var, exception);
        final double[] sTaylor = res.s1dTaylor.sDbl;
        assertEquals(VarDbl.momentum.maxOrder, sTaylor.length);
        assertEquals(sTaylor[0], sin, REPR_DELTA);
        assertEquals(+cos/1, sTaylor[1], REPR_DELTA);
        assertEquals(-sin/2, sTaylor[2], REPR_DELTA);
        assertEquals(-cos/6, sTaylor[3], REPR_DELTA);
        assertEquals(+sin/24, sTaylor[4], REPR_DELTA);
        assertEquals(+cos/120, sTaylor[5], REPR_DELTA);
    }


    @Test
    public void testSin() 
            throws IOException {
        validateSin(0, 0.1, 0, 3e-4, null);
        validateSin(Math.PI/4, 0.1, 6e-8, 2e-4, null);
        validateSin(Math.PI/2, 0.1, 8e-8, 1e-8, null);
        validateSin(Math.PI/4*3, 0.1, 6e-8, 8e-8, null);
        validateSin(Math.PI, 0.1, 1e-23, 2e-7, null);
        validateSin(-Math.PI/4, 0.1, 6e-8, 8e-8, null);
        validateSin(-Math.PI/2, 0.1, 8e-8, 1e-8, null);
        validateSin(-Math.PI/4*3, 0.1, 6e-8, 8e-8, null);
        validateSin(-Math.PI, 0.1, 6e-8, 2e-4, null);
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
}


