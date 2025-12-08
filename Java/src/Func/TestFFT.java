package Func;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import Stats.Histogram;
import Stats.Noise;
import Stats.Stat;
import Type.InitException;
import Type.NotFiniteException;
import Type.NotMonotonicException;
import Type.NotPositiveException;
import Type.NotReliableException;
import Type.NotStableException;
import Type.VarDbl;
import Type.Taylor;



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


class FFT_Signal {
    final SinSource sinSource;
    final FFT fft;

    final int order;
    final int size;
    final int freq;
    final SignalType signalType;
    final VarDbl[] sWave, sFreq;

    FFT_Signal(FFT_Signal signal) {
        this.sinSource = signal.sinSource;
        this.fft = signal.fft;
        this.size = signal.size;
        this.order = signal.order;
        this.freq = signal.freq;
        this.signalType = signal.signalType;
        this.sWave = signal.sWave;
        this.sFreq = signal.sFreq;
    }

    FFT_Signal(int order, SignalType signalType, int freq, SinSource sinSource) 
                    throws InitException, IOException,
                            NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        this.signalType = signalType;
        this.order = order;
        this.freq = freq;

        this.sinSource = sinSource;
        this.fft = new FFT(sinSource);

        this.size = 1 << order;
        int half = size >> 1;
        if (freq >= half) {
            throw new ArithmeticException(String.format("Sin invalid freq %d for order %d", freq, order));
        }
        sWave = new VarDbl[size << 1];
        sFreq = new VarDbl[size << 1];
        if (signalType == SignalType.Sin) {
            for (int i = 0; i < size; ++i) {
                sWave[(i << 1)] = fft.isin.sin(freq *i, order - 1);
                sWave[(i << 1) + 1] = new VarDbl();
                if (i == freq) {
                    sFreq[(i << 1)] = new VarDbl();
                    sFreq[(i << 1) + 1] = new VarDbl(half);
                } else if (i == (size - freq)) {
                    sFreq[(i << 1)] = new VarDbl();
                    sFreq[(i << 1) + 1] = new VarDbl(- half);
                } else {
                    sFreq[(i << 1)] = new VarDbl();
                    sFreq[(i << 1) + 1] = new VarDbl();
                }
            }
        } else if (signalType == SignalType.Cos) {
            for (int i = 0; i < size; ++i) {
                sWave[(i << 1)] = fft.isin.cos(freq *i, order - 1);
                sWave[(i << 1) + 1] = new VarDbl();
                if ((i == freq) || (i == (size - freq))) {
                    sFreq[(i << 1)] = new VarDbl(half);
                    sFreq[(i << 1) + 1] = new VarDbl();
                } else {
                    sFreq[(i << 1)] = new VarDbl();
                    sFreq[(i << 1) + 1] = new VarDbl();
                }
            }
        } else if (signalType == SignalType.Linear) {
            sWave[0] = new VarDbl();
            sWave[1] = new VarDbl();
            sFreq[0] = new VarDbl(((double) size) * (size - 1) / 2);
            sFreq[1] = new VarDbl();
            for (int i = 1; i < size; ++i) {
                sWave[(i << 1)] = new VarDbl(i);
                sWave[(i << 1) + 1] = new VarDbl();
                sFreq[(i << 1)] = new VarDbl(-half);
                sFreq[(i << 1) + 1] = new VarDbl(-half).multiplyInPlace(fft.isin.cos(i, order))
                                        .multiplyInPlace(Taylor.pow(fft.isin.sin(i, order), -1));
            }
        } else {
            throw new IllegalArgumentException(String.format("Unknown signal %s", signalType));
        }
    }
}



class Measure {
    static final int BINDING = 3, DIVIDS = 5;

    final Map<TestType, Stat> sUncertainty = new HashMap<>();
    final Map<TestType, Stat> sValue = new HashMap<>();
    final Map<TestType, Histogram> sHisto = new HashMap<>();

    Measure() {
        for (TestType test: TestType.values()) {
            sUncertainty.put(test, new Stat());
            sValue.put(test, new Stat());
            sHisto.put(test, new Histogram(BINDING, DIVIDS));
        }
    }
}


class FFT_Order extends FFT_Signal {
    static final double NORMALIZED_ERROR_OUTLIER = 1e14;
    static final Map<Integer, Map<NoiseType, Map<Double, Map<SinSource, Measure>>>> ssssAggr = new HashMap<>();

    final Measure measure = new Measure();
    final Measure aggr;

    final NoiseType noiseType;
    final double noise;
    static Noise rand = new Noise();

    final VarDbl[] sFwrd, sBack;
    final VarDbl[] sSpec, sRound, sRev;
    final VarDbl[][] ssSpecStep, ssRoundStep, ssRevStep;

    FFT_Order(FFT_Signal signal, NoiseType noiseType, double noise, boolean traceSteps) 
            throws InitException, IOException, 
                NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        super(signal);
        this.noiseType = noiseType;
        this.noise = Math.abs(noise);

        sFwrd = new VarDbl[size << 1];
        sBack = new VarDbl[size << 1];
        for (int i = 0; i < (size << 1); ++i) { 
            if (noise == 0) {
                sFwrd[i] = new VarDbl(sWave[i]);
                sBack[i] = new VarDbl(sFreq[i]);
            } else {
                sFwrd[i] = new VarDbl(sWave[i]).addInPlace(new VarDbl(getNoise(), noise));
                sBack[i] = new VarDbl(sFreq[i]).addInPlace(new VarDbl(getNoise(), noise));
            }
        }

        sSpec = fft.transform(sFwrd, true, traceSteps);
        if (traceSteps) {
            ssSpecStep = fft.ssStep; 
            ssSpecStep[order + FFT.EXTRA_STEPS - 2] = sFreq;
            ssSpecStep[order + FFT.EXTRA_STEPS - 1] = new VarDbl[size << 1];
        } else {
            ssSpecStep = null;
        }
        sRound = fft.transform(sSpec, false, traceSteps);
        if (traceSteps) {
            ssRoundStep = fft.ssStep; 
            ssRoundStep[order + FFT.EXTRA_STEPS - 2] = sFwrd;
            ssRoundStep[order + FFT.EXTRA_STEPS - 1] = new VarDbl[size << 1];
        } else {
            ssRoundStep = null;
        }
        sRev = fft.transform(sBack, false, traceSteps);
        if (traceSteps) {
            ssRevStep = fft.ssStep; 
            ssRevStep[order + FFT.EXTRA_STEPS - 2] = sWave;
            ssRevStep[order + FFT.EXTRA_STEPS - 1] = new VarDbl[size << 1];
        } else {
            ssRevStep = null;
        } 

        Measure aggr = null;
        if (signalType == SignalType.Sin || signalType == SignalType.Cos) {
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
            aggr = sAggr.get(sinSource);
            if (aggr == null) {
                aggr = new Measure();
                sAggr.put(sinSource, aggr);
            }
        }
        this.aggr = aggr;

        for (int i = 0; i < (size << 1); ++i) {
            accum(TestType.Forward, sSpec[i], sFreq[i], traceSteps? ssSpecStep[order + FFT.EXTRA_STEPS - 1] : null, i);
            accum(TestType.Roundtrip, sRound[i], sFwrd[i], traceSteps? ssRoundStep[order + FFT.EXTRA_STEPS - 1] : null, i);
            accum(TestType.Reverse, sRev[i], sWave[i], traceSteps? ssRevStep[order + FFT.EXTRA_STEPS - 1] : null, i);
        }
    }

    FFT_Order(FFT_Signal signal, NoiseType noiseType, double noise) 
            throws InitException, IOException, 
                NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        this(signal, noiseType, noise, false);
    }

    private double getNoise() {
        switch (noiseType) {
            case Gaussian:
                return rand.gaussian(noise);
            case White:
                return rand.white(noise);
            default:
                return 0;
        }
    }

    private void accum(final TestType test, final VarDbl actual, final VarDbl exptected,  
                       final VarDbl sStep[], int index) 
            throws InitException {
        final VarDbl err = actual.minus(exptected);
        accum(test, actual, err, index, measure);
        if (aggr != null)
            accum(test, actual, err, index, aggr);
        if (sStep != null)
            sStep[index] = err;
    }

    private void accum(final TestType test, final VarDbl actual, final VarDbl error, int index,
            final Measure measure) {
        measure.sUncertainty.get(test).accum(actual.uncertainty(), index);
        measure.sValue.get(test).accum(error.value(), index);
        if (error.uncertainty() > 0) {
            final double norm = error.value()/error.uncertainty();
            if (Math.abs(norm) < NORMALIZED_ERROR_OUTLIER) 
                measure.sHisto.get(test).accum(norm, index);
            else if (sinSource != SinSource.Lib && actual.value() != 0 && actual.uncertainty() != 0)
                System.out.println(String.format("For signal=%s freq=%d noiseType=%s noise=%e test=%s index=%d, normaliized error outlier %e for error %s and actual %s",
                        signalType, freq, noiseType, noise, test, index, norm, error, actual));
        }
    }

    static private void dump(final FileWriter fw, SinSource sinSource, NoiseType noiseType, double noise,
              SignalType signal, int order, int freq, final Measure measure) throws IOException {
        for (TestType test: TestType.values()) {
            fw.write(String.format("%s\t%s\t%.1e\t%s\t%d\t%d\t%s", 
                     sinSource, noiseType, noise, signal, order, freq, test));
            Stat stat = measure.sUncertainty.get(test);
            fw.write(String.format("\t%d\t%.15e\t%.15e\t%.15e\t%d\t%.15e\t%d", 
                    stat.count(), stat.avg(), stat.dev(), stat.min(), stat.minAt(), stat.max(), stat.maxAt()));
            stat = measure.sValue.get(test);
            fw.write(String.format("\t%d\t%.15e\t%.15e\t%.15e\t%d\t%.15e\t%d", 
                    stat.count(), stat.avg(), stat.dev(), stat.min(), stat.minAt(), stat.max(), stat.maxAt()));
            final Histogram histo = measure.sHisto.get(test);
            fw.write(String.format("\t%d\t%.15e\t%.15e\t%.15e\t%d\t%.15e\t%d", 
                    histo.stat().count(), histo.stat().avg(), histo.stat().dev(), histo.stat().min(), histo.stat().minAt(), histo.stat().max(), histo.stat().maxAt()));
            fw.write(String.format("\t%d\t%d", histo.lower(), histo.upper()));
            final double[] sHisto = histo.histo();
            if (sHisto != null) {
                for (int i = 0; i < sHisto.length; ++i) 
                    fw.write(String.format("\t%.6g", sHisto[i]));
            }
            fw.write("\n");
        }
        fw.flush();
    }

    static private void dump(final FileWriter fw, SinSource sinSource, NoiseType noiseType, double noise,
              SignalType signalType, int order, int freq) 
            throws IOException, ArithmeticException, InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        final FFT_Signal signal = new FFT_Signal(order, signalType, freq, sinSource);
        final FFT_Order out = new FFT_Order(signal, noiseType, noise);
        dump(fw, sinSource, noiseType, noise, signalType, order, freq, out.measure);
    }

    static String header() {
        StringBuilder hdrBuilder = new StringBuilder(
                "SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest" +
                "\tUncertainty Count\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Minimum At\tUncertainty Maximum\tUncertainty Maximum At" +
                "\tValue Count\tValue Mean\tValue Deviation\tValue Minimum\tValue Minimum At\tValue Maximum\tValue Maximum At" +
                "\tError Count\tError Mean\tError Deviation\tError Minimum\tError Minimum At\tError Maximum\tError Maximum At" +
                "\tLower\tUpper");
        for (int i = - Measure.BINDING * Measure.DIVIDS; i <= Measure.BINDING * Measure.DIVIDS; ++i) {
            hdrBuilder.append(String.format("\t%.1f", ((double) i) / Measure.DIVIDS));  
        }
        return hdrBuilder.toString();
    }

    static boolean dump(String dumpPath, int minOrder, int maxOrder, List<Integer> sFreq,
            List<SinSource> sSinSource, List<NoiseType> sNoiseType, List<Double> sNoise) {
        final Map<Integer, Map<NoiseType, Map<Double, Map<SinSource, Set<TestType>>>>> sssssAggr = read(dumpPath);
        try (final FileWriter fw = new FileWriter(dumpPath, sssssAggr != null)) {
            if (sssssAggr == null)
                fw.write(FFT_Order.header() + "\n");
            for (NoiseType noiseType: sNoiseType) {
                for (int order = minOrder; order < maxOrder; ++order) {
                    for (final SinSource sinSource : sSinSource) {
                        System.out.println(String.format("%s: Starting noiseType=%s order=%d sinSource=%s", 
                                LocalDateTime.now(), noiseType, order, sinSource));
                        for (final double noise : sNoise) {
                            if (sssssAggr != null) {
                                final Map<NoiseType, Map<Double, Map<SinSource, Set<TestType>>>> ssssAggr = sssssAggr.get(order);
                                if (ssssAggr != null) {
                                    final Map<Double, Map<SinSource, Set<TestType>>> sssAggr = ssssAggr.get(noiseType);
                                    if (sssAggr != null) {
                                        final Map<SinSource, Set<TestType>> ssAggr = sssAggr.get(noise);
                                        if (ssAggr != null) {
                                            final Set<TestType> sAggr = ssAggr.get(sinSource);
                                            if ((sAggr != null) && (sAggr.size() == TestType.values().length))
                                                continue;
                                        }
                                    }
                                }
                            }
                            dump(fw, sinSource, noiseType, noise, SignalType.Linear, order, 0);
                            for (final int freq : sFreq) {
                                if (freq >= (1 << (order - 1)))
                                    continue;
                                dump(fw, sinSource, noiseType, noise, SignalType.Sin, order, freq);
                                dump(fw, sinSource, noiseType, noise, SignalType.Cos, order, freq);
                            }        
                            final Map<NoiseType, Map<Double, Map<SinSource, Measure>>> sssAggr = FFT_Order.ssssAggr.get(order);
                            if (sssAggr == null) {
                                fail(String.format("Found no aggregated result for order=%d", order));
                            }
                            final Map<Double, Map<SinSource, Measure>> ssAggr = sssAggr.get(noiseType);
                            if (ssAggr == null) {
                                fail(String.format("Found no aggregated result for order=%d noiseType=%s", 
                                                    order, noiseType));
                            }
                            final Map<SinSource, Measure> sAggr = ssAggr.get(noise);
                            if (sAggr == null) {
                                fail(String.format("Found no aggregated result for order=%d noiseType=%s noise=%e", 
                                                    order, noiseType, noise));
                            }
                            final Measure aggr = sAggr.get(sinSource);
                            if (aggr == null) {
                                fail(String.format("Found no aggregated result for order=%d noiseType=%s noise=%e sinSource=%s", 
                                order, noiseType, noise, sinSource));
                            }
                            dump(fw, sinSource, noiseType, noise, SignalType.Aggr, order, 0, aggr);
                            fw.flush();
                        }
                    }
                }
            }
        } catch (IOException | ArithmeticException | NotFiniteException | NotReliableException | NotMonotonicException | NotStableException | NotPositiveException | InitException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    static boolean dump(String dumpPath, int maxOrder) {
        return dump(dumpPath, IndexSin.MIN_ORDER + 1, maxOrder, 
            List.of(1, 2,3,4,5,6,7), 
            List.of(SinSource.Prec, SinSource.Quart, SinSource.Lib), 
            List.of(NoiseType.Gaussian, NoiseType.White),
            List.of(0., 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 
                    1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.));
    }

    static boolean dump(String dumpPath) {
        return dump(dumpPath, IndexSin.MAX_ORDER + 1);
    }
    
    static Map<Integer, Map<NoiseType, Map<Double, Map<SinSource, Set<TestType>>>>> read(String dumpPath) {
        final String HEADER = FFT_Order.header();
        Map<Integer, Map<NoiseType, Map<Double, Map<SinSource, Set<TestType>>>>> sssssAggr = null;
        try (final BufferedReader br = new BufferedReader(new FileReader(dumpPath))) {
            Map<NoiseType, Map<Double, Map<SinSource, Set<TestType>>>> ssssAggr = null;
            Map<Double, Map<SinSource, Set<TestType>>> sssAggr = null;
            Map<SinSource, Set<TestType>> ssAggr = null;
            Set<TestType> sAggr = null;
            String line = br.readLine();
            if (line != null && line.equals(HEADER)) {
                for (line = br.readLine(); line != null; line = br.readLine()) {
                    final String[] sWord = line.trim().split("\t");
                    final SinSource sinSource = SinSource.valueOf(sWord[0]);
                    final NoiseType noiseType = NoiseType.valueOf(sWord[1]);
                    final double noise = Double.parseDouble(sWord[2]);
                    final SignalType signalType = SignalType.valueOf(sWord[3]);
                    final int order = Integer.parseInt(sWord[4]);
                    final double freq = Double.parseDouble(sWord[5]);
                    final TestType testType = TestType.valueOf(sWord[6]);
                    if ((order < IndexSin.MIN_ORDER) || (IndexSin.MAX_ORDER < order))
                        throw new IllegalArgumentException(String.format("invalid order %d in line: %s", order, line));
                    if ((noise < 0) || (1 < noise))
                        throw new IllegalArgumentException(String.format("invalid noise %e in line: %s", noise, line));
                    if ((freq < 0) || (7 < freq))
                        throw new IllegalArgumentException(String.format("invalid freq %e in line: %s", freq, line));
                    if (signalType != SignalType.Aggr)
                        continue;
                    if (sssssAggr == null)
                        sssssAggr = new HashMap<>();
                    ssssAggr = sssssAggr.get(order);
                    if (ssssAggr == null) {
                        ssssAggr = new HashMap<>();
                        sssssAggr.put(order, ssssAggr);
                    }
                    sssAggr = ssssAggr.get(noiseType);
                    if (sssAggr == null) {
                        sssAggr = new HashMap<>();
                        ssssAggr.put(noiseType, sssAggr);
                    }
                    ssAggr = sssAggr.get(noise);
                    if (ssAggr == null) {
                        ssAggr = new HashMap<>();
                        sssAggr.put(noise, ssAggr);
                    }
                    sAggr = ssAggr.get(sinSource);
                    if (sAggr == null) {
                        sAggr = new HashSet<>();
                        ssAggr.put(sinSource, sAggr);
                    }
                    sAggr.add(testType);
                }
            } else
                return null;
        } catch (IOException | IllegalArgumentException ex) {
            ex.printStackTrace();
            return null;
        }
        return sssssAggr;
    }
}


class FFT_Step extends FFT_Order {
    FFT_Step(FFT_Signal signal, NoiseType noiseType, double noise) 
            throws InitException, IOException, 
                NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        super(signal, noiseType, noise, true);
    }

    static private void dump(final FileWriter fw, final String context, final VarDbl[] sData) throws IOException {
        for (int imag: List.of(0, 1)) {
            for (int value: List.of(1, 0)) {
                fw.write(String.format("\n%s\t%d\t%d", context, imag, value));
                for (int i = imag; i < sData.length; i += 2) {
                     fw.write(String.format("\t%.17e", (value == 1)? sData[i].value() : sData[i].uncertainty()));
                }
            }
        }
    }

    static private void dump(final FileWriter fw, String context, final TestType testType, 
            final VarDbl[][] ssStep, final VarDbl[] sExpected) 
            throws IOException, InitException {
        for (int step = 0; step < ssStep.length; ++step) {
            if (ssStep[step] == null) {
                final int order = ssStep.length - FFT.EXTRA_STEPS;
                if (step == order + 2)
                    break;
                else
                    throw new RuntimeException(String.format("Invalid termination step %d for %s", step, context));
            }
            final String tag = String.format("%s\t%s\t%d", context, testType, step);
            dump(fw, tag, ssStep[step]);
        }
    }

    static private void dump(final FileWriter fw, SinSource sinSource, int order,
            SignalType signalType, int freq, NoiseType noiseType, double noise) 
            throws IOException, InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        final FFT_Signal signal = new FFT_Signal(order, signalType, freq, sinSource);
        final FFT_Step out = new FFT_Step(signal, noiseType, noise);
        final String context = String.format("%s\t%s\t%.1e\t%s\t%d\t%d", 
                sinSource, noiseType, noise, signalType, order, freq);
        dump(fw, context, TestType.Forward, out.ssSpecStep, out.sFreq);
        dump(fw, context, TestType.Roundtrip, out.ssRoundStep, out.sFwrd);
        dump(fw, context, TestType.Reverse, out.ssRevStep, out.sWave);
    }

    static boolean dump(String dumpPath, SinSource sinSource, int order, 
            NoiseType noiseType, double noise) {
        if (dumpPath == null)
            dumpPath = String.format("FFT_Step_%d_%s.txt", order, sinSource);
        try (final FileWriter fw = new FileWriter(dumpPath)) {
            fw.write("SinSource\tNoiseType\tNoise\tSignal\tOrder\tFreq\tTest\tStep\tImag\tValue");
            final int size = 1 << order;
            final VarDbl[] sCosSin = new VarDbl[size << 1];
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%d", i));
            final IndexSin isin = new IndexSin(sinSource);
            for (int i = 0; i < size; ++i) {
                sCosSin[(i << 1)] = isin.cos(i, order);
                sCosSin[(i << 1) + 1] = isin.sin(i, order);
            }
            final String context = String.format("%s\t\t0\t\t%d\t\tCosSin\t", 
                    sinSource, order);
            dump(fw, context, sCosSin);
            for (int freq = 1; freq < 8 && freq < (size >> 1); ++freq)
                dump(fw, sinSource, order, SignalType.Sin, freq, noiseType, noise);
            for (int freq = 1; freq < 8 && freq < (size >> 1); ++freq)
                dump(fw, sinSource, order, SignalType.Cos, freq, noiseType, noise);
            dump(fw, sinSource, order, SignalType.Linear, 0, noiseType, noise);
       } catch (IOException | InitException | ArithmeticException e) {
            e.printStackTrace();
            return false;
        } catch (NotFiniteException | NotReliableException | NotMonotonicException | NotStableException | NotPositiveException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }
}

public class TestFFT {
    static final int[] sFreq = new int[]{1, 2, 3, 4, 5, 6, 7};
    static final double[] sNoise = new double[]{0, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1};

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
    public void testFFTOrder2Sin() throws InitException, IOException {
        final VarDbl[] sFwrd = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(1),new VarDbl(), 
            new VarDbl(),new VarDbl(), new VarDbl(-1),new VarDbl()};
        final VarDbl[] sSpec = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(2), 
            new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(-2)}; 
        for (SinSource sinSource: SinSource.values()) {
            final FFT fft = new FFT(sinSource);
            validate( sSpec, fft.transform(sFwrd, true), 0);
            validate( sFwrd, fft.transform(sSpec, false), 0);
        }
    }

    @Test
    public void testFFTOrder2Cos() throws InitException, IOException {
        final VarDbl[] sFwrd = new VarDbl[] {
            new VarDbl(1),new VarDbl(), new VarDbl(),new VarDbl(), 
            new VarDbl(-1),new VarDbl(), new VarDbl(0),new VarDbl()};
        final VarDbl[] sSpec = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(2),new VarDbl(), 
            new VarDbl(),new VarDbl(), new VarDbl(2),new VarDbl()}; 
        for (SinSource sinSource: SinSource.values()) {
            final FFT fft = new FFT(sinSource);
            validate( sSpec, fft.transform(sFwrd, true), 0);
            validate( sFwrd, fft.transform(sSpec, false), 0);
        }
    }

    @Test
    public void testFFTOrder3Sin() throws InitException, IOException {
        final double q = Math.sqrt(0.5);
        final VarDbl[] sFwrd = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(q),new VarDbl(),  new VarDbl(1),new VarDbl(), new VarDbl(q),new VarDbl(), 
            new VarDbl(),new VarDbl(), new VarDbl(-q),new VarDbl(), new VarDbl(-1),new VarDbl(), new VarDbl(-q),new VarDbl()};
        final VarDbl[] sSpec = new VarDbl[] {
            new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(4), new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(),
            new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(), new VarDbl(),new VarDbl(-4)}; 
        for (SinSource sinSource: SinSource.values()) {
            final FFT fft = new FFT(sinSource);
            final VarDbl[] sFor = fft.transform(sFwrd, true);
            validate(sSpec, sFor, Math.ulp(2));
            final VarDbl[] sRev = fft.transform(sSpec, false);
            validate(sFwrd, sRev, Math.ulp(2));
        }

        try {
            final FFT_Signal signal = new FFT_Signal(3, SignalType.Sin, 1, SinSource.Lib);
            final FFT_Order out = new FFT_Order(signal, NoiseType.Gaussian, 0);
            assertEquals(1.8679476, out.measure.sHisto.get(TestType.Forward).stat().dev(), 0.01);
        } catch (NotFiniteException | NotReliableException | NotMonotonicException | NotStableException
                | NotPositiveException | InitException | IOException e) {
            e.printStackTrace();
            fail(e.getMessage());
        }

    }

    @Test
    public void testFFTOrder4Sin() throws InitException, IOException {
        final int size = 1 << 4;
        final int freq = 4;
        final double f = Math.PI*2 *freq/size;
        final VarDbl[] sFwrd = new VarDbl[size << 1];
        final VarDbl[] sSpec = new VarDbl[size << 1];
        for (int i = 0; i < size; ++i) {
            sFwrd[(i << 1)] = new VarDbl(Math.sin(f*i));
            sFwrd[(i << 1) + 1] = new VarDbl();
            sSpec[(i << 1)] = new VarDbl();
            if (i == freq)
                sSpec[(i << 1) + 1] = new VarDbl(size/2);
            else if ((size - i) == freq)
                sSpec[(i << 1) + 1] = new VarDbl(-size/2);
            else
                sSpec[(i << 1) + 1] = new VarDbl();
        }
        for (SinSource sinSource: SinSource.values()) {
            final FFT fft = new FFT(sinSource);
            final VarDbl[] sFor = fft.transform(sFwrd, true);
            validate(sSpec, sFor, Math.ulp(2));
            final VarDbl[] sRev = fft.transform(sSpec, false);
            validate(sFwrd, sRev, Math.ulp(2));
        }
    }
  

    @Test 
    public void dump_Step_2_6() {
        for (int order = 2; order < 6; ++order) {
            for (SinSource sinSource: SinSource.values()) {
                final String dumpPath = String.format("./Java/Output/FFT_Step_%d_%s.txt", order, sinSource);
                FFT_Step.dump(dumpPath, sinSource, order, NoiseType.Gaussian, 0);
            }
        }
    }

    @Test 
    public void dump_Order_2_6() {
       FFT_Order.dump("./Java/Output/FFT_Order_2_6.txt", 6);
    }

    // very time comsuming
    @Test 
    public void dump_Order_2_19() throws ArithmeticException, InitException {
        FFT_Order.dump("./Java/Output/FFT_2_19.txt");
    }

    
}
