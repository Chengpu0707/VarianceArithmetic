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

import Func.FFT.SinSource;
import Stats.Histogram;
import Stats.Noise;
import Stats.Stat;
import Type.InitException;
import Type.VarDbl;



/*
 * A class to facilitate testing of FFT
 */
enum SignalType {
    Sin,
    Cos,
    Linear,

    Aggr,  // aggregated Sin and Cos for all frequencies 
} 

enum NoiseType {
    Gaussian,
    White
}

enum TestType {
    Forward,
    Reverse,
    Roundtrip,
}

class Measure {
    static final int BINDING = 3, DIVIDS = 5;

    final Map<TestType, Stat> sStat = new HashMap<>();
    final Map<TestType, Histogram> sHisto = new HashMap<>();

    Measure() {
        for (TestType test: TestType.values()) {
            sStat.put(test, new Stat());
            sHisto.put(test, new Histogram(BINDING, DIVIDS));
        }
    }
}


class Signal {
   static final Map<Integer, Map<NoiseType, Map<Double, Map<SinSource, Measure>>>> ssssAggr = new HashMap<>();

    final SinSource sinType;
    final FFT fft;

    final int order;
    final int size;
    final int freq;
    final SignalType signal;
    final double[] sWave, sFreq;

    final NoiseType noiseType;
    final double noise;
    static Noise rand = new Noise();

    final VarDbl[] sData, sSpec, sRound, sBack, sRev;

    final Measure measure = new Measure();
     Signal(int order, int freq, SignalType signal, NoiseType noiseType, double noise, FFT.SinSource sinType) 
                throws ArithmeticException, InitException, IOException {
        this.order = order;
        this.freq = freq;
        this.signal = signal;
        this.noiseType = noiseType;
        this.noise = Math.abs(noise);

        this.sinType = sinType;
        this.fft = new FFT(sinType);

        size = 1 << order;
        if (freq * 2 >= size) {
            throw new ArithmeticException(String.format("Sin invalid freq %d for order %d", freq, order));
        }
        sWave = new double[size << 1];
        sFreq = new double[size << 1];
        sData = new VarDbl[size << 1];
        sBack = new VarDbl[size << 1];
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
                    sWave[i << 1] = fft.sin(freq * (long) i, order);
                }
                sFreq[(freq << 1) + 1] = peak;
                sFreq[((size - freq) << 1) + 1] = - peak;
                break;
            case Cos:
                for (int i = 0; i < size; ++i) {
                    sWave[i << 1] = fft.cos(freq * (long) i, order);
                }
                sFreq[freq << 1] = peak;
                sFreq[(size - freq) << 1] = peak;
                break;
            case Linear:
                sFreq[0] = size * (size - 1) / 2;
                for (int i = 1; i < size; ++i) {
                    sWave[i << 1] = i;
                    sFreq[i << 1] = -peak;
                    if (order < IndexSin.MAX_ORDER)
                        sFreq[(i << 1) + 1] = -peak / fft.sin(i, order + 1) * fft.cos(i, order + 1);
                    else
                        sFreq[(i << 1) + 1] = -peak / Math.sin(Math.PI * i / size) * Math.cos(Math.PI * i / size);
                }
                break;
             default:
                throw new ArithmeticException(String.format("Unknown signal %s", signal));
        }

        for (int i = 0; i < (size << 1); ++i) { 
            if (noise == 0) {
                sData[i] = new VarDbl(sWave[i] + getNoise());
                sBack[i] = new VarDbl(sFreq[i] + getNoise());
            } else {
                sData[i] = new VarDbl(sWave[i] + getNoise(), noise);
                sBack[i] = new VarDbl(sFreq[i] + getNoise(), noise);
            }
        }

        sSpec = fft.transform(sData, true, 
                        String.format("./Java/Output/testFFTOrder%d%s%d_for.txt", order, signal, freq));
        sRound = fft.transform(sSpec, false,
                        String.format("./Java/Output/testFFTOrder%d%s%d_rnd.txt", order, signal, freq));
        sRev = fft.transform(sBack, false, 
                        String.format("./Java/Output/testFFTOrder%d%s%d_rev.txt", order, signal, freq));

        Measure aggr = null;
        if (signal != SignalType.Linear) {
            Map<NoiseType, Map<Double, Map<SinSource, Measure>>> sssAggr = ssssAggr.get(order);
            if (sssAggr == null) {
                sssAggr = new HashMap<>();
                ssssAggr.put(order, sssAggr);
            }
            Map<Double, Map<SinSource, Measure>> ssAggr = sssAggr.get(noiseType);
            if (ssAggr == null) {
                ssAggr = new HashMap<>();
                sssAggr.put(noiseType, ssAggr);
            }
            Map<SinSource, Measure> sAggr = ssAggr.get(noise);
            if (sAggr == null) {
                sAggr = new HashMap<>();
                ssAggr.put(noise, sAggr);
            }
            aggr = sAggr.get(sinType);
            if (aggr == null) {
                aggr = new Measure();
                sAggr.put(sinType, aggr);
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

    static final FFT fft = new FFT(FFT.SinSource.IndexSin);

    @Test
    public void testBitReverseIndices() {
        assertArrayEquals( new int[] {0, 2, 1, 3}, FFT.bitReversedIndices(2));
        assertArrayEquals( new int[] {0, 4, 2, 6, 1, 5, 3, 7}, FFT.bitReversedIndices(3));
        assertArrayEquals( new int[] {0, 8, 4, 12, 2, 10, 6, 14, 
                                      1, 9, 5, 13, 3, 11, 7, 15}, FFT.bitReversedIndices(4));
                                      
        for (int order = 2; order <= IndexSin.MAX_ORDER; ++order) {
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

    void validate(final VarDbl[] sExpected, final VarDbl[] sRes, double delta) {
        assertEquals( sExpected.length, sRes.length);
        for (int i = 0; i < sExpected.length; ++i)
            assertEquals(sExpected[i].value(), sRes[i].value(), Math.ulp(2));
    }

    @Test
    public void testFFTOrder2Sin() {
        VarDbl[] sData = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(1),new VarDbl(), 
            new VarDbl(),new VarDbl(), new VarDbl(-1),new VarDbl()};
        VarDbl[] sSpec = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(2), 
            new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(-2)}; 
        try {
            validate( sSpec, fft.transform(sData, true, "./Java/Output/testFFTOrder2Sin_for.txt"), 0);
            validate( sData, fft.transform(sSpec, false, "./Java/Output/testFFTOrder2Sin_rev.txt"), 0);
        } catch (InitException | IOException e) {
            fail(e.getMessage());
        }       
    }

    @Test
    public void testFFTOrder2Cos() {
        VarDbl[] sData = new VarDbl[] {
            new VarDbl(1),new VarDbl(), new VarDbl(),new VarDbl(), 
            new VarDbl(-1),new VarDbl(), new VarDbl(0),new VarDbl()};
        VarDbl[] sSpec = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(2),new VarDbl(), 
            new VarDbl(),new VarDbl(), new VarDbl(2),new VarDbl()}; 
        try {
            validate( sSpec, fft.transform(sData, true, "./Java/Output/testFFTOrder2Sin_for.txt"), 0);
            validate( sData, fft.transform(sSpec, false, "./Java/Output/testFFTOrder2Sin_rev.txt"), 0);
        } catch (InitException | IOException e) {
            fail(e.getMessage());
        }       
    }

    @Test
    public void testFFTOrder3Sin() {
        final double q = Math.sqrt(0.5);
        try {
            final VarDbl[] sData = new VarDbl[] {
                new VarDbl(),new VarDbl(), new VarDbl(q),new VarDbl(),  new VarDbl(1),new VarDbl(), new VarDbl(q),new VarDbl(), 
                new VarDbl(),new VarDbl(), new VarDbl(-q),new VarDbl(), new VarDbl(-1),new VarDbl(), new VarDbl(-q),new VarDbl()};
            final VarDbl[] sSpec = new VarDbl[] {
                new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(4), new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(),
                new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(-4)}; 
            final VarDbl[] sFor = fft.transform(sData, true, "./Java/Output/testFFTOrder3Sin_for.txt");
            validate(sSpec, sFor, Math.ulp(2));
            final VarDbl[] sRev = fft.transform(sSpec, false, "./Java/Output/testFFTOrder3Sin_rev.txt");
            validate(sData, sRev, Math.ulp(2));
        } catch (InitException | IOException e) {
            fail(e.getMessage());
        }       
    }

    @Test
    public void testFFTOrder4Sin() {
        final int size = 1 << 4;
        final int freq = 4;
        final double f = Math.PI*2 *freq/size;
        try {
            final VarDbl[] sData = new VarDbl[size << 1];
            final VarDbl[] sSpec = new VarDbl[size << 1];
            for (int i = 0; i < size; ++i) {
                sData[(i << 1)] = new VarDbl(Math.sin(f*i));
                sData[(i << 1) + 1] = new VarDbl();
                sSpec[(i << 1)] = new VarDbl();
                if (i == freq)
                    sSpec[(i << 1) + 1] = new VarDbl(size/2);
                else if ((size - i) == freq)
                    sSpec[(i << 1) + 1] = new VarDbl(-size/2);
                else
                    sSpec[(i << 1) + 1] = new VarDbl();
            }
            final VarDbl[] sFor = fft.transform(sData, true, "./Java/Output/testFFTOrder4Sin_for.txt");
            validate(sSpec, sFor, Math.ulp(2));
            final VarDbl[] sRev = fft.transform(sSpec, false, "./Java/Output/testFFTOrder4Sin_rev.txt");
            validate(sData, sRev, Math.ulp(2));

            final Signal signal = new Signal(4, freq, SignalType.Sin, NoiseType.White, 0, FFT.SinSource.IndexSin);
            validate(sSpec, signal.sSpec, Math.ulp(2));
            validate(sData, signal.sRev, Math.ulp(2));
            validate(sData, signal.sRound, Math.ulp(2));
            final Measure measure = Signal.ssssAggr.get(4).get(NoiseType.White).get(0.).get(FFT.SinSource.IndexSin);
            assertEquals(0, measure.sHisto.get(TestType.Forward).stat().count());
        } catch (InitException | IOException e) {
            fail(e.getMessage());
        }

    }
    
    static void dump(final FileWriter fw, final String tag, TestType test, String part, 
                     VarDbl[] sOut, double[] sIdeal) throws IOException {
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

    static void inspect(double sigma) {
        final String pathOut = "./Java/Output/FFTVar_spec.txt";
        final int maxOrder = 7;
        try (final FileWriter fw = new FileWriter(pathOut)) {
            fw.write("SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest\tI/O\tMeasure\tPart");
            for (int i = 0; i < (2 << maxOrder); ++i)
                fw.write(String.format("\t%d", i));
            for (FFT.SinSource sinType: FFT.SinSource.values()) {
                for (final NoiseType noiseType: NoiseType.values()) {
                    for (int n = 0; n < sNoise.length; ++n) {
                        final  double noise = sNoise[n];
                        for (final SignalType signal: SignalType.values()) {
                            if (signal == SignalType.Aggr) 
                                continue;
                            for (int order = 4; order <= maxOrder; ++order) {
                                for (int f = 0; f < sFreq.length; ++f) {
                                    final int freq = sFreq[f];
                                    Signal out = new Signal(order, freq, signal, noiseType, noise, sinType);
                                    final String tag = String.format("%s\t%s\t%.3e\t%s\t%d\t%d", 
                                            sinType, noiseType, noise, signal, order, freq);
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
            }
        } catch (InitException e) {
            fail(e.getMessage());
        } catch (ArithmeticException e) {
            fail(e.getMessage());
        } catch (IOException e) {
            fail(e.getMessage());
        }
    }

    // TODO: @Test
    public void inspect() {
        inspect(0);
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

    void dump(final FileWriter fw, SinSource sinType, NoiseType noiseType, double noise,
              SignalType signal, int order, int freq, final Measure measure) throws IOException {
        for (TestType test: TestType.values()) {
            fw.write(String.format("%s\t%s\t%.1e\t%s\t%d\t%d\t%s", 
                     sinType, noiseType, noise, signal, order, freq, test));
            dump(fw, measure.sStat.get(test), measure.sHisto.get(test));
        }
    }

    void run(final FileWriter fw, SinSource sinType, NoiseType noiseType, double noise,
              SignalType signal, int order, int freq) 
            throws IOException, ArithmeticException, InitException {
        final Signal out = new Signal(order, freq, signal, noiseType, noise, sinType);
        dump(fw, sinType, noiseType, noise, signal, order, freq, out.measure);
    }

 
    //TODO: generate FFT @Test
    public void dump() {
        final String pathOut = String.format("./Java/Output/FFT_%d_%d.txt", 4, IndexSin.MAX_ORDER);
        try (final FileWriter fw = new FileWriter(pathOut)) {
            fw.write("SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest");
            fw.write("\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Maximum");
            fw.write("\tError Mean\tError Deviation\tError Minimum\tError Maximum");
            for (int i = - Measure.BINDING * Measure.DIVIDS; i <= Measure.BINDING * Measure.DIVIDS; ++i) {
                fw.append(String.format("\t%.1f", ((double) i) / Measure.DIVIDS));  
            }
            fw.write("\n");
            for (int order = 4; order <= IndexSin.MAX_ORDER; ++order) {
                for (final double noise : sNoise) {
                    for (final NoiseType noiseType: NoiseType.values()) {
                        for (final SinSource sinType : SinSource.values()) {
                            run(fw, sinType, noiseType, noise, SignalType.Linear, order, 0);
                            for (final int freq : sFreq) {
                                run(fw, sinType, noiseType, noise, SignalType.Sin, order, freq);
                                run(fw, sinType, noiseType, noise, SignalType.Cos, order, freq);
                            }        
                        }
                    }
                }
                final Map<NoiseType, Map<Double, Map<SinSource, Measure>>> sssAggr = Signal.ssssAggr.get(order);
                if (sssAggr == null) {
                    fail(String.format("Found no aggregated result for order=%d", order));
                }
                for (NoiseType noiseType: NoiseType.values()) {
                    final Map<Double, Map<SinSource, Measure>> ssAggr = sssAggr.get(noiseType);
                    if (ssAggr == null) {
                        fail(String.format("Found no aggregated result for order=%d noiseType=%s", 
                                            order, noiseType));
                    }
                    for (double noise: sNoise) {
                        final Map<SinSource, Measure> sAggr = ssAggr.get(noise);
                        if (sAggr == null) {
                            fail(String.format("Found no aggregated result for order=%d noiseType=%s noise=%e", 
                                                order, noiseType, noise));
                        }
                        for (final SinSource sinType: SinSource.values()) {
                            final Measure aggr = sAggr.get(sinType);
                            if (aggr == null) {
                                fail(String.format("Found no aggregated result for order=%d noiseType=%s noise=%e sinType=%s", 
                                order, noiseType, noise, sinType));
                            }
                            dump(fw, sinType, noiseType, noise, SignalType.Aggr, order, 0, aggr);
                        }
                    }
                }
            }
        } catch (ArithmeticException | InitException | IOException e) {
            fail(e.getMessage());
        }
    }

}
