package Func;


import Type.InitException;
import Type.VarDbl;


public class FFT {
    /*
     * A cache of the bit-reversed index
     */
    static final int EXTRA_STEPS = 5;
        // 0: Input
        // 1: bit reverse
        // 2 to order + 1: intermediate steps
        // order + 2: output
        // order + 3: exptected
        // order + 4: error

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

    public final IndexSin isin;
    VarDbl[][] ssStep;

    public FFT(SinSource sinSource) {
        this.isin = new IndexSin(sinSource);;
    }

    void addStep(VarDbl[] sRes, int step) {
        ssStep[step] = new VarDbl[sRes.length];
        for (int i = 0; i < sRes.length; ++i)
            ssStep[step][i] = new VarDbl(sRes[i]);
    }

    /*
     * 1-dimentional Fast Fourier Transformation (FFT)
     * 
     * @param sInput    an array of size (2<<order), with each datum contains (real, image)
     * @param forward   true for forware transformation, false for backward transformation.
     * 
     * @return      an array of size (2<<order), with each datum contains (real, image)
     */
    public VarDbl[] transform(final VarDbl[] sInput, boolean forward, boolean traceSteps) 
            throws InitException {
        int order = IndexSin.MIN_ORDER;
        for (; order <= IndexSin.MAX_ORDER; ++order) {
            if ((2 << order) == sInput.length) {
                break;
            }
        }
        if (order > IndexSin.MAX_ORDER) {
            throw new IllegalArgumentException(
                String.format("Invalid input array size=%d for fft.transform()", sInput.length));
        }

        final VarDbl[] sRes = new VarDbl[sInput.length];
        if (traceSteps) {
            ssStep = new VarDbl[order + EXTRA_STEPS][];
            ssStep[0] = sInput;
        }
       
        final int[] sIndex = bitReversedIndices(order);
        for (int i = 0; i < sIndex.length; i++) {
            final int j = sIndex[i];
            sRes[(i << 1)] = sInput[j << 1].clone();
            sRes[(i << 1) + 1] = sInput[(j << 1) + 1].clone();
        }
        if (traceSteps)
            addStep(sRes, 1);

        for (int i = 0; i < (sIndex.length - 1); i += 2 ) {
            final VarDbl rt = sRes[(i << 1)].clone(), it = sRes[(i << 1) + 1].clone();
            sRes[(i << 1)].addInPlace(sRes[(i << 1) + 2]);
            sRes[(i << 1) + 1].addInPlace(sRes[(i << 1) + 3]);
            sRes[(i << 1) + 2] = rt.minusInPlace(sRes[(i << 1) + 2]);
            sRes[(i << 1) + 3] = it.minusInPlace(sRes[(i << 1) + 3]);
        }
        if (traceSteps)
            addStep(sRes, 2);

        for (int o = 1, k = 4; o < order; ++o, k <<= 1) {
            for (int j = 0; j < (k >> 1); j++) {
                final VarDbl cos = isin.cos(j, o);
                final VarDbl sin = isin.sin(forward? j : -j, o);
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
            if (traceSteps)
                addStep(sRes, o + 2);
        }   // for o
        
        if (!forward) {
            for (int i = 0; i < sRes.length; ++i  ) {
                sRes[i].multiplyInPlace(new VarDbl(1.0/(1L << order), 0));
            }
        }
        if (traceSteps)
            addStep(sRes, order + 2);

        return sRes;
    }

    public VarDbl[] transform(final VarDbl[] sInput, boolean forward) 
                throws InitException {
        return transform(sInput, forward, false);
    }

}

