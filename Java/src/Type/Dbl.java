package Type;

import Type.IReal.ValueException;



/*
 * A tool class to decompose double into different components:
 *    Dbl(double): to decompose a double into (exp, val, neg), with true value as (neg? -1 : +1) * val * 2^exp
 *    toDouble(): to compose a double from (exp, val, neg)
 * Dbl is expected to hold finite value only, so both Dbl(double) and toDouble() may through Type.IReal.ValueException.
 * 
 * exp is stored as a signed value with an offset of Dbl.DOUBLE_EXP_OFFSET
 * 
 * The significand has an extra bit Dbl.DOUBLE_VAL_EXTRA when the stored exp is not (Dbl.DOUBLE_EXP_MIN - 1).
 */
public class Dbl {
    int exp;   // exponent, in signed offset format
    long val;    // significand
    boolean neg; // sign

    static final long DOUBLE_SIGN_MASK = 1L << (Double.SIZE - 1);
    static final int  DOUBLE_EXP_SHIFT = 52;
    static final long DOUBLE_EXP_MASK = ((1L << (Double.SIZE - 1 - Dbl.DOUBLE_EXP_SHIFT)) - 1) << Dbl.DOUBLE_EXP_SHIFT;
    static final long DOUBLE_EXP_OFFSET = Double.MAX_EXPONENT + Dbl.DOUBLE_EXP_SHIFT;
    static final int  DOUBLE_EXP_MAX = Double.MAX_EXPONENT - Dbl.DOUBLE_EXP_SHIFT;
    static final int  DOUBLE_EXP_MIN = (int) - DOUBLE_EXP_OFFSET + 1;
    static final long DOUBLE_VAL_EXTRA = 1L << Dbl.DOUBLE_EXP_SHIFT;
    static final long DOUBLE_VAL_MASK = Dbl.DOUBLE_VAL_EXTRA - 1;
    static final long DOUBLE_VAL_MAX = Dbl.DOUBLE_VAL_EXTRA + Dbl.DOUBLE_VAL_MASK;

    public Dbl( final int exp, final boolean neg, final long val ) throws ValueException {
        if (val < 0) {
            throw new IReal.ValueException(String.format("Init Dbl with negative val %d", val));
        }
        this.neg = neg;
        this.exp = exp;
        this.val = val;
    }
    public Dbl( final Dbl dbl ) {
        this.neg = dbl.neg;
        this.exp = dbl.exp;
        this.val = dbl.val;
    }

    public Dbl(final double value) throws ValueException {
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException("Init Dbl with NaN");
        }
        final long bits = Double.doubleToLongBits(value);
        neg = (bits & Dbl.DOUBLE_SIGN_MASK) == Dbl.DOUBLE_SIGN_MASK;
        exp = (short) (((bits & Dbl.DOUBLE_EXP_MASK) >>> Dbl.DOUBLE_EXP_SHIFT) - Dbl.DOUBLE_EXP_OFFSET);
        val = bits & Dbl.DOUBLE_VAL_MASK;
        if (exp == Dbl.DOUBLE_EXP_MIN - 1) {
            ++exp;
        } else {
            val |= Dbl.DOUBLE_VAL_EXTRA;
        }
    }

    public double toDouble() throws ValueException {
        int exp = this.exp;
        long val = this.val;
        if (exp < Dbl.DOUBLE_EXP_MIN) {
            final Round rnd = new Round(val);
            while ((rnd.val > 0) && (exp < Dbl.DOUBLE_EXP_MIN)) {
                rnd.upOnce();
                ++exp;
            }
            val = rnd.val;
        }
        if (val == 0) {
            return 0;
        } 
        if (val > Dbl.DOUBLE_VAL_MAX) {
            final Round rnd = new Round(val);
            while (rnd.val > Dbl.DOUBLE_VAL_MAX) {
                rnd.upOnce();
                ++exp;
            }
            val = rnd.val;
        } else {
            while ((val < DOUBLE_VAL_EXTRA) && (exp > (Dbl.DOUBLE_EXP_MIN))) {
                val <<= 1;
                --exp;
            }
            if (val < DOUBLE_VAL_EXTRA) {
                --exp;
            }
        }
        if (exp > Dbl.DOUBLE_EXP_MAX) {
            throw new ValueException(String.format("val * 2^(%d) is not finite for Dbl", this.val, this.exp));
        }
        return Double.longBitsToDouble( (neg? Dbl.DOUBLE_SIGN_MASK : 0L) | 
            ((exp + Dbl.DOUBLE_EXP_OFFSET) << DOUBLE_EXP_SHIFT) | 
            (val & Dbl.DOUBLE_VAL_MASK));
    }
}
