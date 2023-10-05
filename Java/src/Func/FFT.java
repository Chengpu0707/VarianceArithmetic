package Func;

public class FFT {
    /*
     * The maximual order of FFT calculation.
     * The size of FFT is limited to 2^MAX_ORDER
     */
    private static final int MAX_ORDER = 16;

    /*
    * A cached sine function with resolution down to PI/2^(MAX_ORDER+1)
    */
    private static final int MAX_SIZE = (1 << MAX_ORDER);
    private static final double[] sSin = new double[MAX_SIZE + 1];
    static {
        for (int i = 0; i <= MAX_SIZE; ++i) {
            sSin[i] = Math.sin(Math.PI * 0.5 * i / MAX_SIZE);
        }
    }

    static double sin(int idx, int order) {
        if (order < 2) {
             throw new IllegalArgumentException(
                String.format("The order %d < 2", order, MAX_ORDER));
        }
        if (MAX_ORDER + 2 < order) {
            throw new IllegalArgumentException(
                String.format("The order %d is too large for the cached max order %d", order, MAX_ORDER));
        }
        if (idx < 0) {
            return -sin(-idx, order);
        }
        int u = idx << (MAX_ORDER - (order - 2));
        final int quart = u >> MAX_ORDER;
        u -= quart << MAX_ORDER;
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
                    String.format("The internal error for geting sine for index %d and order%d/%d", idx, order, MAX_ORDER));
        }
    }

    static double cos(int idx, int order) {
        if (order < 2) {
             throw new IllegalArgumentException(
                String.format("The order %d < 2", order, MAX_ORDER));
        }
        if (MAX_ORDER + 2 < order) {
            throw new IllegalArgumentException(
                String.format("The order %d is too large for the cached max order %d", order, MAX_ORDER));
        }
        if (idx < 0) {
            return -cos(-idx, order);
        }
        int u = idx << (MAX_ORDER - (order - 2));
        final int quart = u >> MAX_ORDER;
        u -= quart << MAX_ORDER;
        switch (quart % 4) {
            case 0:
                return sSin[MAX_SIZE - u];
            case 1: 
                return -sSin[u];
            case 2:
                return -sSin[MAX_SIZE - u];
            case 3:
                return sSin[u];
            default:
                throw new UnknownError(
                    String.format("The internal error for geting sine for index %d and order%d/%d", idx, order, MAX_ORDER));
        }
    }
}
