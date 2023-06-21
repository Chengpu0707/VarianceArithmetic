package Type;

import Type.IReal.ValueException;



/*
 * A tool class to decompose double into different components:
 *    Dbl(double): to decompose a double into (exp, val, neg), with true value as (neg? -1 : +1) * val * 2^exp
 *    toDouble(): to compose a double from (exp, val, neg)
 * Dbl is expected to hold finite value only, so both Dbl(double) and toDouble() may through Type.IReal.ValueException.
 *  *) exp is stored as a signed value (with offset of Dbl.DOUBLE_EXP_OFFSET already subtracted)
 *  *) The significand has an extra bit Dbl.DOUBLE_VAL_EXTRA when the stored exp is not (Dbl.DOUBLE_EXP_MIN - 1).
 * 
 * Dbl also has a member rndErr, to minimize rounding error, with the following functions:
 *      upOnce()
 *      toExp(exp)
 */
public class Dbl {
    int exp;        // exponent
    long val;       // significand
    boolean neg;    // sign
    boolean rndErr; // round error for val

    static final long DOUBLE_SIGN_MASK = 1L << (Double.SIZE - 1);
    static final int  DOUBLE_EXP_SHIFT = 52;
    static final long DOUBLE_EXP_MASK = ((1L << (Double.SIZE - 1 - Dbl.DOUBLE_EXP_SHIFT)) - 1) << Dbl.DOUBLE_EXP_SHIFT;
    static final long DOUBLE_EXP_OFFSET = Double.MAX_EXPONENT + Dbl.DOUBLE_EXP_SHIFT;
    static final int  DOUBLE_EXP_MAX = Double.MAX_EXPONENT - Dbl.DOUBLE_EXP_SHIFT;
    static final int  DOUBLE_EXP_MIN = (int) - DOUBLE_EXP_OFFSET + 1;
    static final long DOUBLE_VAL_EXTRA = 1L << Dbl.DOUBLE_EXP_SHIFT;
    static final long DOUBLE_VAL_MASK = Dbl.DOUBLE_VAL_EXTRA - 1;
    static final long DOUBLE_VAL_MAX = Dbl.DOUBLE_VAL_EXTRA + Dbl.DOUBLE_VAL_MASK;

    public Dbl( final int exp, final boolean neg, final long val, final boolean rndErr ) throws ValueException {
        if (val < 0) {
            throw new IReal.ValueException(String.format("Init Dbl with negative val %d", val));
        }
        this.neg = neg;
        this.exp = exp;
        this.val = val;
        this.rndErr = rndErr;
    }
    public Dbl( final Dbl other ) {
        this.neg = other.neg;
        this.exp = other.exp;
        this.val = other.val;
        this.rndErr = other.rndErr;
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

    static final long VAL_BITS = 62;
    static final long VAL_MAX =  (1L << VAL_BITS) - 1;
    static final long VAL_EXTRA = 1L << (VAL_BITS - 1);
    static final long[] BYTES = new long[7];
    static {
        for (int i = 0; i < 7; ++i) {
            BYTES[i] = 1L << (8 * i);
        }
    }

    /*
     * round up once
     */
    void upOnce() {
        if ((val & 1) != 0) {
            if (rndErr) {
                val += 1;
                rndErr = false;
            } else {
                rndErr = true;
            }
        }
        val >>= 1;
        ++exp;
    }

    /*
     * Round up from this.exp to larger exp.  There is no limitation.
     * 
     * Round down from this.exp to smaller exp.  The exp is limited to keep val < VAL_MAX.
     * 
     * Return ture if exp is reached.
     */
    boolean toExp(int exp) {
        int shift = exp - this.exp;
        if (shift <= 0) {
            if (val == 0) {
                this.exp = exp;
                return true;
            }
            if (val >= VAL_EXTRA) {
                return false;
            }
            shift = -shift;
            int i = 1;
            for (; (i < 7) && (val > BYTES[i]); ++i) {
                continue;
            }
            if (i < 7) {
                if (shift < (VAL_BITS - 8 * i)) {
                    val <<= shift;
                    this.exp = exp;
                    return true;
                }
                shift -= (VAL_BITS - 8 * i);
                val <<= (VAL_BITS - 8 * i);
                this.exp -= (VAL_BITS - 8 * i);
            }
            for (; (shift > 0) && (val < VAL_EXTRA); --shift) {
                val <<= 1;
                --this.exp;
            }
            return (this.exp == exp);
        }

        this.exp = exp;
        final long msb = 1L << shift;
        final long remain = val & (msb - 1);
        final long half = msb >> 1;
        val >>>= shift;
        if ((remain > half) || ((remain == half) && rndErr)) {
            val += 1;
            rndErr = false;
        } else if (remain > 0) {
            rndErr = true;
        }
        return true;
    }
}
