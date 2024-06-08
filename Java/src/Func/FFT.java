package Func;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import Func.FFT.SinSource;
import Stats.Histogram;
import Stats.Noise;
import Stats.Stat;
import Type.InitException;
import Type.VarDbl;

public class FFT {

    enum SinSource {
        IndexSin,
        LibSin
    }

    private static final IndexSin isin = new IndexSin(IndexSin.MAX_ORDER);

    /*
     * A cache of the bit-reversed index
     */
    static int[][] sBitReversedIndex = new int[IndexSin.MAX_ORDER - 1][];
    static {
        for (int order = 2; order <= IndexSin.MAX_ORDER; ++order) {
            final int[] sRes = sBitReversedIndex[order - 2] = new int[1 << order];
            final int N = 1 << order;
            final int M = N >> 1;
            int i, j, k;
            j = 0;
            for (i = 0; i < N; i++) {
                sRes[i] = j;
                // next j from NumericalRecipesinC.pdf
                k = M; 
                while ((k != 0) && (j >= k)) {
                    j -= k;
                    k >>= 1;
                }
                j += k;
            }
        }
    }

    /*
     * Generate the bit reversed index of "order": from NumericalRecipesinC.pdf
     */
    static int[] bitReversedIndices(int order) {
        if ((order < 2) || ((IndexSin.MAX_ORDER) < order))
            return null;
        return sBitReversedIndex[order - 2];
    }

    private final SinSource _sinType;

    public FFT(SinSource sinType) {
        _sinType = sinType;
    }

    public double sin(long freq, int order) {
        switch (_sinType) {
        case IndexSin:
            if (order < IndexSin.MIN_ORDER) {
                throw new IllegalArgumentException(
                   String.format("The order %d < %d for fft.sin()", order, IndexSin.MIN_ORDER));
           }
           if (isin.order() < order) {
               throw new IllegalArgumentException(
                   String.format("The order %d > %d for fft.sin()", order, isin.order()));
           }
           return isin.sin(freq *(1L <<(IndexSin.MAX_ORDER - order + 1)));
        case LibSin:            
            return Math.sin(Math.PI *freq /(1L << (order - 1)));
        default:
            throw new IllegalArgumentException(
                String.format("Unknown SinSource %s for fft.sin()", _sinType));
        }
    }

    public double cos(long freq, int order) {
        switch (_sinType) {
        case IndexSin:
            if (order < IndexSin.MIN_ORDER) {
                throw new IllegalArgumentException(
                   String.format("The order %d < %d for fft.cos()", order, IndexSin.MIN_ORDER));
           }
           if (isin.order() < order) {
               throw new IllegalArgumentException(
                   String.format("The order %d > %d for fft.cos()", order, isin.order()));
           }
           return isin.cos(freq *(1L <<(IndexSin.MAX_ORDER - order + 1)));
        case LibSin:            
            return Math.cos(Math.PI *freq /(1L << (order - 1)));
        default:
            throw new IllegalArgumentException(
                String.format("Unknown SinSource %s for fft.cos()", _sinType));
        }
    }

    /*
     * 1-dimentional Fast Fourier Transformation (FFT)
     * 
     * @param sData     an array of size (2<<order), with each datum contains (real, image)
     * @param forward   true for forware transformation, false for backward transformation.
     * 
     * @return      an array of size (2<<order), with each datum contains (real, image)
     */
    public VarDbl[] transform(final VarDbl[] sData, boolean forward, final String dumpPath) 
            throws InitException, IOException {
        int order = IndexSin.MIN_ORDER;
        for (; order <= IndexSin.MAX_ORDER; ++order) {
            if ((2 << order) == sData.length) {
                break;
            }
        }
        if (order > IndexSin.MAX_ORDER) {
            throw new IllegalArgumentException(
                String.format("Invalid input array size %d for fft.transform()", sData.length));
        }
        final int size = 1 << order;

        FileWriter fw = (dumpPath == null)? null : new FileWriter(dumpPath);
        if (fw != null) {
            fw.write("Step\tPart");
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%d Value\t%d Uncertainty", i, i));
            fw.write(String.format("\nInput\tReal"));
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%.6e\t%.3e", 
                        sData[(i << 1)].value(), sData[(i << 1)].uncertainty()));
            fw.write("\nInput\tImag");
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%.6e\t%.3e", 
                        sData[(i << 1) + 1].value(), sData[(i << 1) + 1].uncertainty()));
            fw.flush(); 
        }

        final VarDbl[] sRes = new VarDbl[2 << order];
        
        final int[] sIndex = bitReversedIndices(order);
        for (int i = 0; i < sIndex.length; i++) {
            final int j = sIndex[i];
            sRes[(i << 1)] = sData[j << 1].clone();
            sRes[(i << 1) + 1] = sData[(j << 1) + 1].clone();
        }
        if (fw != null) {
            fw.write("\n0\tReal");
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%.6e\t%.3e", 
                        sRes[(i << 1)].value(), sRes[(i << 1)].uncertainty()));
            fw.write("\n0\tImag");
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%.6e\t%.3e", 
                        sRes[(i << 1) + 1].value(), sRes[(i << 1) + 1].uncertainty()));
            fw.flush(); 
        }

        for (int i = 0; i < (sIndex.length - 1); i += 2 ) {
            final VarDbl rt = sRes[(i << 1)].clone(), it = sRes[(i << 1) + 1].clone();
            sRes[(i << 1)].add(sRes[(i << 1) + 2], true);
            sRes[(i << 1) + 1].add(sRes[(i << 1) + 3], true);
            sRes[(i << 1) + 2] = rt.minus(sRes[(i << 1) + 2], true);
            sRes[(i << 1) + 3] = it.minus(sRes[(i << 1) + 3], true);
        }
        if (fw != null) {
            fw.write("\n1\tReal");
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%.6e\t%.3e", 
                        sRes[(i << 1)].value(), sRes[(i << 1)].uncertainty()));
            fw.write("\n1\tImag");
            for (int i = 0; i < size; ++i)
                fw.write(String.format("\t%.6e\t%.3e", 
                        sRes[(i << 1) + 1].value(), sRes[(i << 1) + 1].uncertainty()));
            fw.flush(); 
        }

        for (int o = 2, k = 4; o <= order; ++o, k <<= 1) {
            for (int j = 0; j < (k >> 1); j++) {
                final VarDbl cos = new VarDbl(cos( j, o ));
                final VarDbl sin = new VarDbl(forward? sin( j, o ) : -sin( j, o ));
                for (int i = 0; i < sIndex.length; i += k ) {
                    final int idx0 = (i + j) << 1;
                    final int idx1 = idx0 + k;
                    final VarDbl r1 = sRes[idx1].clone();
                    final VarDbl i1 = sRes[idx1 + 1].clone();
            
                    final VarDbl rd = r1.multiply(cos).minus(i1.multiply(sin));
                    final VarDbl id = i1.multiply(cos).add(r1.multiply(sin));
            
                    sRes[idx1] = sRes[idx0].minus(rd);
                    sRes[idx1 + 1] = sRes[idx0 + 1].minus(id);
                    sRes[idx0].add(rd, true);
                    sRes[idx0 + 1].add(id, true);
                }   // for( i
            }
            if (fw != null) {
                fw.write(String.format("\n%d\tReal", o));
                for (int i = 0; i < size; ++i)
                    fw.write(String.format("\t%.6e\t%.3e", 
                            sRes[(i << 1)].value(), sRes[(i << 1)].uncertainty()));
                fw.write(String.format("\n%d\tImag", o));
                for (int i = 0; i < size; ++i)
                    fw.write(String.format("\t%.6e\t%.3e", 
                            sRes[(i << 1) + 1].value(), sRes[(i << 1) + 1].uncertainty()));
                fw.flush(); 
            }
        }
        
        if (!forward) {
            for (int i = 0; i < (sIndex.length << 1); i ++ ) {
                sRes[i].multiply(new VarDbl(1.0/(1L << order), 0), true);
            }
            if (fw != null) {
                fw.write("\nOutput\tReal");
                for (int i = 0; i < size; ++i)
                    fw.write(String.format("\t%.6e\t%.3e", 
                            sRes[(i << 1)].value(), sRes[(i << 1)].uncertainty()));
                fw.write("\nOutput\tImag");
                for (int i = 0; i < size; ++i)
                    fw.write(String.format("\t%.6e\t%.3e", 
                            sRes[(i << 1) + 1].value(), sRes[(i << 1) + 1].uncertainty()));
                fw.flush(); 
            }
        }
        if (fw != null)
            fw.close();

        return sRes;
    }
    public VarDbl[] transform(final VarDbl[] sData, boolean forward) 
            throws InitException {
        try {
            return transform(sData, forward, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }
    public VarDbl[] transform(final double[] sData, boolean forward, final String dumpPath) 
            throws InitException, IOException {
        final VarDbl[] sInput = new VarDbl[sData.length];
        for (int i = 0; i < sData.length; ++i)
            sInput[i] = new VarDbl(sData[i]);
        return transform(sInput, forward, dumpPath);
    }
    public VarDbl[] transform(final double[] sData, boolean forward) 
            throws InitException {
        try {
            return transform(sData, forward, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }
}

/*
 * A class to facilitate testing of FFT
 */
class Signal {

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

    static class Measure {
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
    static final Map<Integer, Map<NoiseType, Map<Double, Map<SinSource, Measure>>>> ssssAggr = new HashMap<>();

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




