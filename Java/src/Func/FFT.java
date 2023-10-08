package Func;

import Type.IReal;

public class FFT {
    /*
     * The maximual order of FFT calculation.
     * The size of FFT is limited to 2^MAX_QUARD_ORDER
     */
    private static final int MAX_QUARD_ORDER = 16;

    /*
    * A cached sine function with resolution down to PI/2^(MAX_QUARD_ORDER+1)
    */
    private static final int MAX_SIZE = (1 << MAX_QUARD_ORDER);
    private static final double[] sSin = new double[MAX_SIZE + 1];
    static {
        for (int i = 0; i <= MAX_SIZE; ++i) {
            sSin[i] = Math.sin(Math.PI * 0.5 * i / MAX_SIZE);
        }
    }

    /*
     * @return: sine(idx/(1<<order)*2*PI)
     * @param order: the size of [0, 2*PI] is (1<<order), so 2 < order <= MAX_QUARD_ORDER+2
     * @param idx: the index in [0, 2*PI]
     */
    static double sin(int idx, int order) {
        if (order < 2) {
             throw new IllegalArgumentException(
                String.format("The order %d < 2", order, MAX_QUARD_ORDER));
        }
        if (MAX_QUARD_ORDER + 2 < order) {
            throw new IllegalArgumentException(
                String.format("The order %d is too large for the cached max order %d", order, MAX_QUARD_ORDER));
        }
        if (idx < 0) {
            return -sin(-idx, order);
        }
        int u = idx << (MAX_QUARD_ORDER - (order - 2));
        final int quart = u >> MAX_QUARD_ORDER;
        u -= quart << MAX_QUARD_ORDER;
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
                    String.format("The internal error for geting sine for index %d and order%d/%d", idx, order, MAX_QUARD_ORDER));
        }
    }

    static double cos(int idx, int order) {
        if (order < 2) {
             throw new IllegalArgumentException(
                String.format("The order %d < 2", order, MAX_QUARD_ORDER));
        }
        return sin(idx + (1 << (order - 2)), order);
    }

    /*
     * @param sData: an array of size (2<<order), with each datum contains (real, image)
     */
    static IReal[] transform(final IReal[] sData, boolean forward) {
        int order = 2;
        for (; order <= MAX_QUARD_ORDER + 2; ++order) {
            if ((2 << order) == sData.length) {
                break;
            }
        }
        if (order > MAX_QUARD_ORDER + 2) {
            return null;
        }
        final IReal[] sRes = new IReal[2 << order];
        
        return null;
    }


}
