package Func;

/*
 * A cached sine function with resolution down to PI /2^order().
 * It has integer input to sin(), cos() and tan() to avoid rounding errors.
 */
public class IndexSin {
    /*
     * The maximual order of FFT calculation.
     * The size of FFT is limited to 2^MAX_QUARD_ORDER
     */
    static public final int MIN_ORDER = 2;
    static public final int MAX_ORDER = 18;

    private final int _order;
    private final int _half;
    private final double[] _sSin;

    public IndexSin(int order) {
        if (order < MIN_ORDER) {
             throw new IllegalArgumentException(
                String.format("The order %d < %d for IndexSin", order, MIN_ORDER));
        }
        if (MAX_ORDER < order) {
            throw new IllegalArgumentException(
                String.format("The order %d > %d for IndexSin", order, MAX_ORDER));
        }
        _order = order;
        _half = 1 << (order - 1);
        _sSin = new double[_half + 1];
        final int size = _half / 2;
        for (int i = 0; i < size; ++i) {
            _sSin[i] = Math.sin(Math.PI*0.5 * i/_half);
            _sSin[i + size] = Math.cos(Math.PI*0.5 * (size - i)/_half);
        }
        _sSin[_half] = 1;
    }

    public int order() {
        return _order;
    }

    int get_index(long freq) {
        boolean neg = freq < 0;
        long div = Math.abs(freq) / this._half;
        long rem = Math.abs(freq) % _half;
        if ((div & 1) == 1) {
            rem = _half - rem;
        }
        if ((div & 2) == 2) {
            neg = !neg;
        }
        return (int) (neg? -rem : rem);
    }
    
    public double sin(long freq) {
        final int idx = this.get_index(freq);
        return (idx >= 0)? _sSin[idx] : -_sSin[-idx];
    }

    public double cos(long freq) {
        return sin(freq + _half);
    }
        
    public double tran(long freq) {
        return sin(freq) / cos(freq);
    }
}
