package Func;

import java.io.FileWriter;
import java.io.IOException;

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
            sRes[(i << 1)].addInPlace(sRes[(i << 1) + 2]);
            sRes[(i << 1) + 1].addInPlace(sRes[(i << 1) + 3]);
            sRes[(i << 1) + 2] = rt.minusInPlace(sRes[(i << 1) + 2]);
            sRes[(i << 1) + 3] = it.minusInPlace(sRes[(i << 1) + 3]);
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
                    sRes[idx0].addInPlace(rd);
                    sRes[idx0 + 1].addInPlace(id);
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
                sRes[i].multiplyInPlace(new VarDbl(1.0/(1L << order), 0));
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

