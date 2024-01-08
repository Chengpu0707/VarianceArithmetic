package Func;

import Stats.Noise;
import Type.IReal;
import Type.IntvDbl;
import Type.VarDbl;
import Type.IReal.TypeException;
import Type.IReal.UncertaintyException;
import Type.IReal.ValueException;

public class FFT {
    /*
     * The maximual order of FFT calculation.
     * The size of FFT is limited to 2^MAX_QUARD_ORDER
     */
    private static final int MAX_QUARD_ORDER = 16;
    static final int MAX_ORDER = MAX_QUARD_ORDER + 2;

    /*
    * A cached sine function with resolution down to PI/2^(MAX_QUARD_ORDER+1)
    * Use "noiseDev" to reduce systemic numerical errors in sin()
    */
    private static final int MAX_SIZE = (1 << MAX_QUARD_ORDER);
    private static final double[] sSin = new double[MAX_SIZE + 1];
    private static final double noiseDev = 0;
    static {
        Noise noise = new Noise();
        for (int i = 0; i <= MAX_SIZE; ++i) {
            sSin[i] = Math.sin(Math.PI * 0.5 * i / MAX_SIZE) + noise.white(noiseDev);
        }
    }

    /*
     * @return: sine(idx/(1<<order)*2*PI)
     * @param order: the size of [0, 2*PI] is (1<<order), so 2 < order <= MAX_QUARD_ORDER+2
     * @param idx: the index in [0, 2*PI]
     */
    static double sin(long idx, int order) {
        if (idx < 0) {
            return -sin(-idx, order);
        }
        if (order < 2) {
             throw new IllegalArgumentException(
                String.format("The order %d < 2", order, MAX_QUARD_ORDER));
        }
        if (MAX_QUARD_ORDER + 2 < order) {
            throw new IllegalArgumentException(
                String.format("The order %d is too large for the cached max order %d", order, MAX_QUARD_ORDER));
        }
        final int size = 1 << (order - 2);
        int u = ((int) (idx % size)) << (MAX_ORDER - order);
        final int quart = (int) (idx / size);
        switch (quart % 4) {
            case 0:
                return sSin[u];
            case 1: 
                return sSin[MAX_SIZE - u];
            case 2:
                return -sSin[u];
            case 3:
                return -sSin[MAX_SIZE - u];
            default:
                throw new UnknownError(
                    String.format("Invalid quart %d for geting sine for index %d and order %d/%d", 
                                  quart, idx, order, MAX_QUARD_ORDER));
        }
    }

    static double cos(long idx, int order) {
        if (order < 2) {
             throw new IllegalArgumentException(
                String.format("The order %d < 2", order, MAX_QUARD_ORDER));
        }
        return sin(idx + (1 << (order - 2)), order);
    }

    /*
     * A cache of the bit-reversed index
     */
    static int[][] sBitReversedIndex = new int[MAX_QUARD_ORDER + 1][];
    static {
        for (int order = 2; order <= MAX_QUARD_ORDER + 2; ++order) {
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
        if ((order < 2) || ((MAX_QUARD_ORDER + 2) < order))
            return null;
        return sBitReversedIndex[order - 2];
    }

    /*
     * @param sData: an array of size (2<<order), with each datum contains (real, image)
     * @param forward: true for forware transformation, false for backward transformation.
     * @param useOriginalArray: change sData without allocate a result array
     */
    static public IReal[] transform(final IReal[] sData, boolean forward, 
                boolean useOriginalArray) 
            throws TypeException, ValueException, UncertaintyException {
        int order = 2;
        for (; order <= MAX_QUARD_ORDER + 2; ++order) {
            if ((2 << order) == sData.length) {
                break;
            }
        }
        if (order > MAX_ORDER) {
            throw new IllegalArgumentException(String.format("Invalid input array size %d", sData.length));
        }
        final RealType realType;
        if (sData[0] instanceof VarDbl) 
            realType = RealType.Var;
        else if (sData[0] instanceof IntvDbl)
            realType = RealType.Intv;
        else
            throw new IllegalArgumentException(String.format("Invalid input array type %s", sData[0]));

        final IReal[] sRes = useOriginalArray? sData : new IReal[2 << order];
        
        final int[] sIndex = bitReversedIndices(order);
        for (int i = 0; i < sIndex.length; i++) {
            final int j = sIndex[i];
            if (useOriginalArray) {
                if (i > j) {
                    IReal x = sRes[(i << 1)];
                    sRes[(i << 1)] = sRes[(j << 1)];
                    sRes[(j << 1)] = x;
                    x = sRes[(i << 1) + 1];
                    sRes[(i << 1) + 1] = sRes[(j << 1) + 1];
                    sRes[(j << 1) + 1] = x;
                }
            } else {
                sRes[(i << 1)] = sData[j << 1].clone();
                sRes[(i << 1) + 1] = sData[(j << 1) + 1].clone();
            }
        }

        for (int i = 0; i < (sIndex.length - 1); i += 2 ) {
            final IReal rt = sRes[(i << 1)].clone(), it = sRes[(i << 1) + 1].clone();
            sRes[(i << 1)].add(sRes[(i << 1) + 2]);
            sRes[(i << 1) + 1].add(sRes[(i << 1) + 3]);
            sRes[(i << 1) + 2] = rt.minus(sRes[(i << 1) + 2]);
            sRes[(i << 1) + 3] = it.minus(sRes[(i << 1) + 3]);
        }

        for (int o = 2, k = 2; o <= order; ++o, k <<= 1) {
            for (int j = 0; j < k; j++) {
                final IReal cos = (realType == RealType.Var)
                                    ? new VarDbl(cos( j, o ))
                                    : new IntvDbl(cos( j, o ));
                final IReal sin = (realType == RealType.Var)
                                    ? new VarDbl(forward? sin( j, o ) : -sin( j, o ))
                                    : new IntvDbl(forward? sin( j, o ) : -sin( j, o ));
                for (int i = 0; i < sIndex.length; i += (k << 1) ) {
                    final int idx0 = (i + j) << 1;
                    final int idx1 = idx0 + (k << 1);
                    final IReal r1 = sRes[idx1];
                    final IReal i1 = sRes[idx1 + 1];
            
                    final IReal rd = r1.clone().multiply(cos).minus(i1.clone().multiply(sin));
                    final IReal id = i1.clone().multiply(cos).add(r1.clone().multiply(sin));
            
                    sRes[idx1] = sRes[idx0].clone().minus(rd);
                    sRes[idx1 + 1] = sRes[idx0 + 1].clone().minus(id);
                    sRes[idx0].add(rd);
                    sRes[idx0 + 1].add(id);
                }   // for( i
            }
        }
        if (!forward) {
            for (int i = 0; i < (sIndex.length << 1); i ++ ) {
                sRes[i].shift(-order);
            }
        }

        return sRes;
    }

    static public IReal[] transform(final IReal[] sData, boolean forward)
            throws TypeException, ValueException, UncertaintyException {
        return transform( sData, forward, false);
    }
}
