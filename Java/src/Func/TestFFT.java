package Func;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import Stats.Histogram;
import Stats.Noise;
import Stats.Stat;
import Type.IReal;
import Type.VarDbl;
import Type.IntvDbl;
import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;


enum SignalType {
    Sin,
    Cos,
    Linear,

    Aggr,  // aggregated Sin and Cos for all frequencies 
}

class Signal {
    static Noise rand = new Noise();

    final RealType realType;

    final int order;
    final int size;
    final int freq;
    final SignalType signal;
    final double[] sWave, sFreq;

    final NoiseType noiseType;
    final double noise;

    final IReal[] sData, sSpec, sRound, sBack, sRev;

    final Measure measure = new Measure();
    static final Map<RealType, Map<Integer, Map<NoiseType, Map<Double, Measure>>>> ssssAggr = new HashMap<>();

    Signal(RealType realType, int order, int freq, SignalType signal, 
                NoiseType noiseType, double noise) 
            throws ArithmeticException, TypeException, ValueException, UncertaintyException {
        this.realType = realType;
        this.order = order;
        this.freq = freq;
        this.signal = signal;
        this.noiseType = noiseType;
        this.noise = Math.abs(noise);

        size = 1 << order;
        if (freq * 2 >= size) {
            throw new ArithmeticException(String.format("Sin invalid freq %d for order %d", freq, order));
        }
        sWave = new double[size << 1];
        sFreq = new double[size << 1];

        switch (realType) {
            case Var:
                sData = new VarDbl[size << 1];
                sBack = new VarDbl[size << 1];
                break;
        
            case Intv:
                sData = new IntvDbl[size << 1];
                sBack = new IntvDbl[size << 1];
                break;
            default:
                throw new ArithmeticException(String.format("Unknown realType %s", realType));
        }

        switch (noiseType) {
            case Gaussian:
                break;
            case White:
                break;
            default:
                throw new ArithmeticException(String.format("Unknown noise type %s", noiseType));
        }

        final int peak = size >> 1;
        switch (signal) {
            case Sin:
                for (int i = 0; i < size; ++i) {
                    sWave[i << 1] = FFT.sin(freq * (long) i, order);
                }
                sFreq[(freq << 1) + 1] = peak;
                sFreq[((size - freq) << 1) + 1] = - peak;
                break;
            case Cos:
                for (int i = 0; i < size; ++i) {
                    sWave[i << 1] = FFT.cos(freq * (long) i, order);
                }
                sFreq[freq << 1] = peak;
                sFreq[(size - freq) << 1] = peak;
                break;
            case Linear:
                sFreq[0] = size * (size - 1) / 2;
                for (int i = 1; i < size; ++i) {
                    sWave[i << 1] = i;
                    sFreq[i << 1] = -peak;
                    sFreq[(i << 1) + 1] = -peak / Math.tan(Math.PI * i / size);
                }
                break;
             default:
                throw new ArithmeticException(String.format("Unknown signal %s", signal));
        }

        switch (realType) {
            case Var:
                for (int i = 0; i < (size << 1); ++i) {
                    if (noise == 0) {
                        sData[i] = new VarDbl(sWave[i] + getNoise());
                        sBack[i] = new VarDbl(sFreq[i] + getNoise());
                    } else {
                        sData[i] = new VarDbl(sWave[i] + getNoise(), noise);
                        sBack[i] = new VarDbl(sFreq[i] + getNoise(), noise);
                    }
                }
                break;       
            case Intv:
                for (int i = 0; i < (size << 1); ++i) {
                    if (noise == 0) {
                        sData[i] = new IntvDbl(sWave[i] + getNoise());
                        sBack[i] = new IntvDbl(sFreq[i] + getNoise());
                    } else {
                        sData[i] = new IntvDbl(sWave[i] + getNoise(), noise);
                        sBack[i] = new IntvDbl(sFreq[i] + getNoise(), noise);
                    }
                }
                break;
        }

        sSpec = FFT.transform(sData, true);
        sRound = FFT.transform(sSpec, false);
        sRev = FFT.transform(sBack, false);

        Measure aggr = null;
        if (signal != SignalType.Linear) {
            Map<Integer, Map<NoiseType, Map<Double, Measure>>> sssAggr = ssssAggr.get(realType);
            if (sssAggr == null) {
                sssAggr =  new HashMap<>();
                ssssAggr.put(realType, sssAggr);
            }
            Map<NoiseType, Map<Double, Measure>> ssAggr = sssAggr.get(order);
            if (ssAggr == null) {
                ssAggr = new HashMap<>();
                sssAggr.put(order, ssAggr);
            }
            Map<Double, Measure> sAggr = ssAggr.get(noiseType);
            if (sAggr == null) {
                sAggr = new HashMap<>();
                ssAggr.put(noiseType, sAggr);
            }
            aggr = sAggr.get(noise);
            if (aggr == null) {
                aggr = new Measure();
                sAggr.put(noise, aggr);
            }
        }

        for (int i = 0; i < (size << 1); ++i) {
            final double unc1 = sSpec[i].uncertainty();
            measure.sStat.get(TestType.Forward).accum(unc1);
            if (aggr != null)
                aggr.sStat.get(TestType.Forward).accum(unc1);
            if (unc1 > 0) {
                 measure.sHisto.get(TestType.Forward).accum((sSpec[i].value() - sFreq[i])/unc1);
                if (aggr != null)
                    aggr.sHisto.get(TestType.Forward).accum((sSpec[i].value() - sFreq[i])/unc1);
            }

            final double unc2 = sRound[i].uncertainty();
            measure.sStat.get(TestType.Roundtrip).accum(unc2);
            if (aggr != null)
                aggr.sStat.get(TestType.Roundtrip).accum(unc2);
            if (unc2 > 0) {
                measure.sHisto.get(TestType.Roundtrip).accum((sRound[i].value() - sData[i].value())/unc2);
                if (aggr != null)
                    aggr.sHisto.get(TestType.Roundtrip).accum((sRound[i].value() - sData[i].value())/unc2);
            }

            final double unc3 = sRev[i].uncertainty();
            measure.sStat.get(TestType.Reverse).accum(unc3);
            if (aggr != null)
                aggr.sStat.get(TestType.Reverse).accum(unc3);
            if (unc3 > 0) {
                 measure.sHisto.get(TestType.Reverse).accum((sRev[i].value() - sWave[i])/unc3);
                if (aggr != null)
                    aggr.sHisto.get(TestType.Reverse).accum((sRev[i].value() - sWave[i])/unc3);
            }
        }
    }

    double getNoise() {
        switch (noiseType) {
            case Gaussian:
                return rand.gaussian(noise);
            case White:
                return rand.white(noise);
            default:
                return 0;
        }
    }
}




public class TestFFT {
    static final int[] sFreq = new int[]{1, 2, 3, 4, 5, 6, 7};
    static final double[] sNoise = new double[]{0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 
                                                1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1};

    static final double q1 = Math.sin(Math.PI/8);
    static final double q2 = Math.sin(Math.PI/4);
    static final double q3 = Math.sin(Math.PI/8*3);

    @Test
    public void testSine() {
        assertEquals(-q3, FFT.sin(-5, 4), Math.ulp(q3));

        assertEquals(-1,  FFT.sin(-4, 4), Math.ulp(1.0));
        assertEquals(-q3, FFT.sin(-3, 4), Math.ulp(q3));
        assertEquals(-q2, FFT.sin(-2, 4), Math.ulp(q2));
        assertEquals(-q1, FFT.sin(-1, 4), Math.ulp(q1));

        assertEquals(0, FFT.sin(0, 4), Double.MIN_VALUE);
        assertEquals(q1, FFT.sin(1, 4), Math.ulp(q1));
        assertEquals(q2, FFT.sin(2, 4), Math.ulp(q2));
        assertEquals(q3, FFT.sin(3, 4), Math.ulp(q3));

        assertEquals(1, FFT.sin(4, 4), Math.ulp(1.0));
        assertEquals(q3, FFT.sin(5, 4), Math.ulp(q3));
        assertEquals(q2, FFT.sin(6, 4), Math.ulp(q2));
        assertEquals(q1, FFT.sin(7, 4), Math.ulp(q1));

        assertEquals(0, FFT.sin(8, 4), Double.MIN_VALUE);
        assertEquals(-q1, FFT.sin(9, 4), Math.ulp(q1));
        assertEquals(-q2, FFT.sin(10, 4), Math.ulp(q2));
        assertEquals(-q3, FFT.sin(11, 4), Math.ulp(q3));
        assertEquals(-1, FFT.sin(12, 4), Math.ulp(1.0));

        assertEquals(-q3, FFT.sin(13, 4), Math.ulp(q3));
        assertEquals(-q2, FFT.sin(14, 4), Math.ulp(q2));
        assertEquals(-q1, FFT.sin(15, 4), Math.ulp(q1));
        assertEquals(0, FFT.sin(16, 4), Double.MIN_VALUE);

        assertEquals(q1, FFT.sin(17, 4), Math.ulp(q1));
        assertEquals(q2, FFT.sin(18, 4), Math.ulp(q2));
        assertEquals(q3, FFT.sin(19, 4), Math.ulp(q3));
        
        assertEquals(1, FFT.sin(20, 4), Math.ulp(1.0));
    }

    @Test
    public void testLargeIndexSin() {
        assertEquals(FFT.sin(16, FFT.MAX_ORDER), FFT.sin(268451838, 15), Math.ulp(1.0));
        assertEquals(0, FFT.sin(-2147483648L, 17), Math.ulp(1.0));
        assertEquals(FFT.sin(49152L, 18), FFT.sin(2147532800L, 18), Math.ulp(1.0));
    }

    @Test
    public void testCosine() {
        assertEquals(-q1, FFT.cos(-5, 4), Math.ulp(q1));

        assertEquals(0, FFT.cos(-4, 4), Double.MIN_VALUE);
        assertEquals(q1, FFT.cos(-3, 4), Math.ulp(q1));
        assertEquals(q2, FFT.cos(-2, 4), Math.ulp(q2));
        assertEquals(q3, FFT.cos(-1, 4), Math.ulp(q3));

        assertEquals(1, FFT.cos(0, 4), Math.ulp(1.0));
        assertEquals(q3, FFT.cos(1, 4), Math.ulp(q3));
        assertEquals(q2, FFT.cos(2, 4), Math.ulp(q2));
        assertEquals(q1, FFT.cos(3, 4), Math.ulp(q1));

        assertEquals(0, FFT.cos(4, 4), Double.MIN_VALUE);
        assertEquals(-q1, FFT.cos(5, 4), Math.ulp(q1));
        assertEquals(-q2, FFT.cos(6, 4), Math.ulp(q2));
        assertEquals(-q3, FFT.cos(7, 4), Math.ulp(q3));

        assertEquals(-1, FFT.cos(8, 4), Math.ulp(1.0));
        assertEquals(-q3, FFT.cos(9, 4), Math.ulp(q3));
        assertEquals(-q2, FFT.cos(10, 4), Math.ulp(q2));
        assertEquals(-q1, FFT.cos(11, 4), Math.ulp(q1));

        assertEquals(0, FFT.cos(12, 4), Double.MIN_VALUE);
        assertEquals(q1, FFT.cos(13, 4), Math.ulp(q1));
        assertEquals(q2, FFT.cos(14, 4), Math.ulp(q2));
        assertEquals(q3, FFT.cos(15, 4), Math.ulp(q3));

        assertEquals(1, FFT.cos(16, 4), Math.ulp(1.0));
        assertEquals(q3, FFT.cos(17, 4), Math.ulp(q3));
        assertEquals(q2, FFT.cos(18, 4), Math.ulp(q2));
        assertEquals(q1, FFT.cos(19, 4), Math.ulp(q1));

        assertEquals(0, FFT.cos(20, 4), Double.MIN_VALUE);
    }

    @Test
    public void testBitReverseIndices() {
        assertArrayEquals( new int[] {0, 2, 1, 3}, FFT.bitReversedIndices(2));
        assertArrayEquals( new int[] {0, 4, 2, 6, 1, 5, 3, 7}, FFT.bitReversedIndices(3));
        assertArrayEquals( new int[] {0, 8, 4, 12, 2, 10, 6, 14, 
                                      1, 9, 5, 13, 3, 11, 7, 15}, FFT.bitReversedIndices(4));
                                      
        for (int order = 2; order <= FFT.MAX_ORDER; ++order) {
            final int[] sRes = FFT.bitReversedIndices(order);
            assertNotNull(String.format("order=%d", order), sRes);
            for (int i = 0; i < sRes.length; ++i) {
                final int j = sRes[i];
                int br = 0, org = i;
                for (int k = 0; k < order; ++k) {
                    br <<= 1;
                    br |= (org & 1);
                    org >>= 1;
                }
                assertEquals(String.format("order=%d, index=%d: %d != %d", order, i, br, j), br, j);
            }
        }
    }

    static void dump(final FileWriter fw, final String tag, TestType test, String part, 
                     IReal[] sOut, double[] sIdeal) throws IOException {
        fw.write(String.format("\n%s\t%s\t%s\tValue\tReal", tag, test, part));
        for (int i = 0; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", sOut[i].value()));
        fw.write(String.format("\n%s\t%s\t%s\tError\tReal", tag, test, part));
        for (int i = 0; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", Math.abs(sOut[i].value() - sIdeal[i])));
        fw.write(String.format("\n%s\t%s\t%s\tUncertainty\tReal", tag, test, part));
        for (int i = 0; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", sOut[i].uncertainty()));
        fw.write(String.format("\n%s\t%s\t%s\tNormalized\tReal", tag, test, part));
        for (int i = 0; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", (sOut[i].value() == sIdeal[i])?
                        0 : Math.abs(sOut[i].value() - sIdeal[i])/ sOut[i].uncertainty()));

        fw.write(String.format("\n%s\t%s\t%s\tValue\tImag", tag, test, part));
        for (int i = 1; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", sOut[i].value()));
        fw.write(String.format("\n%s\t%s\t%s\tError\tImag", tag, test, part));
        for (int i = 1; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", Math.abs(sOut[i].value() - sIdeal[i])));
        fw.write(String.format("\n%s\t%s\t%s\tUncertainty\tImag", tag, test, part));
        for (int i = 1; i < sOut.length; i += 2)
            fw.write(String.format("\t%e", sOut[i].uncertainty()));
        fw.write(String.format("\n%s\t%s\t%s\tNormalized\tImag", tag, test, part));
        for (int i = 1; i < sOut.length; i += 2) 
            fw.write(String.format("\t%e", (sOut[i].value() == sIdeal[i])?
                        0 : Math.abs(sOut[i].value() - sIdeal[i])/ sOut[i].uncertainty()));
    }

    static void inspect(RealType realType, double sigma) {
        final String pathOut = String.format("./Java/Output/FFT%s_spec.txt", realType);
        final int maxOrder = 7;
        try (final FileWriter fw = new FileWriter(pathOut)) {
            fw.write("RealType\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest\tI/O\tMeasure\tPart");
            for (int i = 0; i < (2 << maxOrder); ++i)
                fw.write(String.format("\t%d", i));
            for (final NoiseType noiseType: NoiseType.values()) {
                for (int n = 0; n < sNoise.length; ++n) {
                    final  double noise = sNoise[n];
                    for (final SignalType signal: SignalType.values()) {
                        if (signal == SignalType.Aggr) 
                            continue;
                        for (int order = 4; order <= maxOrder; ++order) {
                            for (int f = 0; f < sFreq.length; ++f) {
                                final int freq = sFreq[f];
                                Signal out = new Signal(realType, order, freq, signal, noiseType, noise);
                                final String tag = String.format("%s\t%s\t%.3e\t%s\t%d\t%d", 
                                        realType, noiseType, noise, signal, order, freq);
                                dump(fw, tag, TestType.Forward, "input", out.sData, out.sWave);
                                dump(fw, tag, TestType.Forward, "output", out.sSpec, out.sFreq);
                                dump(fw, tag, TestType.Reverse, "input", out.sBack, out.sFreq);
                                dump(fw, tag, TestType.Reverse, "output", out.sRev, out.sWave);
                                
                                if (sigma > 0) {
                                    StringBuilder sb = new  StringBuilder();
                                    for (int i = 0; i < (out.size << 1); ++i) {
                                        if (Math.abs(out.sSpec[i].value() - out.sFreq[i]) > sigma * out.sSpec[i].uncertainty())
                                            sb.append(String.format("%s Forward %d: %.3e vs %.3e\n", tag, i, 
                                                out.sSpec[i].value() - out.sFreq[i], sigma * out.sSpec[i].uncertainty()));
                                        if (Math.abs(out.sRev[i].value() - out.sWave[i]) > sigma * out.sRev[i].uncertainty())
                                            sb.append(String.format("%s Reverse %d: %.3e vs %.3e\n", tag, i, 
                                                out.sRev[i].value() - out.sWave[i], sigma * out.sRev[i].uncertainty()));
                                    }
                                    assertEquals("", sb.toString());
                                }
                            }
                        }
                    }
                }
            }
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        } catch (ArithmeticException e) {
            fail(e.getMessage());
        } catch (IOException e) {
            fail(e.getMessage());
        }
    }

    void inspect(RealType realType) {
        inspect(realType, 0);
    }

    @Test
    public void inspectVarDbl() {
        inspect(RealType.Var);
    }



    void dump(final FileWriter fw, final Stat stat, final Histogram histo) throws IOException {
        fw.write(String.format("\t%.3e\t%.3e\t%.3e\t%.3e", 
                 stat.avg(), stat.dev(), stat.min(), stat.max()));
        fw.write(String.format("\t%.3e\t%.3e\t%.3e\t%.3e", 
                 histo.stat().avg(), histo.stat().dev(), histo.stat().min(), histo.stat().max()));
        final double[] sHisto = histo.histo();
        if (sHisto != null) {
            for (int i = 0; i < sHisto.length; ++i) 
                fw.write(String.format("\t%.3g", sHisto[i]));
        }
        fw.write("\n");
    }

    void dump(final FileWriter fw, RealType realType, NoiseType noiseType, double noise,
              SignalType signal, int order, int freq, final Measure measure) throws IOException {
        for (TestType test: TestType.values()) {
            fw.write(String.format("%s\t%s\t%.1e\t%s\t%d\t%d\t%s", 
                    realType, noiseType, noise, signal, order, freq, test));
            dump(fw, measure.sStat.get(test), measure.sHisto.get(test));
        }
    }

    void run(final FileWriter fw, RealType realType, NoiseType noiseType, double noise,
              SignalType signal, int order, int freq) 
            throws IOException, ArithmeticException, TypeException, ValueException, UncertaintyException {
        final Signal out = new Signal(realType, order, freq, signal, noiseType, noise);
        dump(fw, realType, noiseType, noise, signal, order, freq, out.measure);
    }

    void run(final FileWriter fw, RealType realType, double noise,
              SignalType signal, int order, int freq) 
            throws IOException, ArithmeticException, TypeException, ValueException, UncertaintyException {
        run(fw, realType, NoiseType.Gaussian, noise, SignalType.Sin, order, freq);
        run(fw, realType, NoiseType.Gaussian, noise, SignalType.Cos, order, freq);
        run(fw, realType, NoiseType.White,  noise, SignalType.Sin, order, freq);
        run(fw, realType, NoiseType.White,  noise, SignalType.Cos, order, freq);
    }


    void dump(RealType realType) {
        final String pathOut = String.format("./Java/Output/FFT%s.txt", realType);
        try (final FileWriter fw = new FileWriter(pathOut)) {
            fw.write("RealType\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest");
            fw.write("\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Maximum");
            fw.write("\tError Mean\tError Deviation\tError Minimum\tError Maximum");
            for (int i = - Measure.BINDING * Measure.DIVIDS; i <= Measure.BINDING * Measure.DIVIDS; ++i) {
                fw.append(String.format("\t%.1f", ((double) i) / Measure.DIVIDS));  
            }
            fw.write("\n");
            for (int order = 4; order <= FFT.MAX_ORDER; ++order) {
                final int maxFreq = (1 << (order - 1)) - 1;
                for (final double noise : sNoise) {
                    run(fw, realType, NoiseType.Gaussian, noise, SignalType.Linear, order, 0);
                    run(fw, realType, NoiseType.White, noise, SignalType.Linear, order, 0);
                    run(fw, realType, noise, SignalType.Linear, order, maxFreq);
                }
                for (final int freq : sFreq) {
                    for (final double noise : sNoise) {
                        run(fw, realType, noise, SignalType.Sin, order, freq);
                    }
                }
                final Map<Integer, Map<NoiseType, Map<Double, Measure>>> sssAggr = Signal.ssssAggr.get(realType);
                if (sssAggr == null) {
                    fail(String.format("Found no aggregated result for %s order=%d", 
                         realType, order));
                }
                final Map<NoiseType, Map<Double, Measure>> ssAggr = sssAggr.get(order);
                if (ssAggr == null) {
                    fail(String.format("Found no aggregated result for %s order=%d", 
                         realType, order));
                }
                for (NoiseType noiseType: NoiseType.values()) {
                    final Map<Double, Measure> sAggr = ssAggr.get(noiseType);
                    if (sAggr == null) {
                        fail(String.format("Found no aggregated result for %s order=%d %s", 
                            realType, order, noiseType));
                    }
                    for (double noise: sNoise) {
                        final Measure aggr = sAggr.get(noise);
                        if (aggr == null) {
                            fail(String.format("Found no aggregated result for %s order=%d %s noise=%.3e", 
                                realType, order, noiseType, noise));
                        }
                        dump(fw, realType, noiseType, noise, SignalType.Aggr, order, 0, aggr);
                    }
                }
            }
        } catch (TypeException | ValueException | UncertaintyException e) {
            fail(e.getMessage());
        } catch (ArithmeticException e) {
            fail(e.getMessage());
        } catch (IOException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void dumpVarDbl() {
        dump(RealType.Var);
    }
}
