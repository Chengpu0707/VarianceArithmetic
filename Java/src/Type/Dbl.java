package Type;

import Type.IReal.ValueException;



/*
 * A tool class to decompose double into different components:
 *    Dbl(double): to decompose a double into (exp, val, neg), with true value as (neg? -1 : +1) * val * 2^exp
 *    normalize(): to fit (exp, val, neg) into a double format.
 *    toDouble(): to compose a double from (exp, val, neg).  It also normalizes the object
 * 
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
        final long msb = Double.doubleToLongBits(value);
        neg = (msb & Dbl.DOUBLE_SIGN_MASK) == Dbl.DOUBLE_SIGN_MASK;
        exp = (short) (((msb & Dbl.DOUBLE_EXP_MASK) >>> Dbl.DOUBLE_EXP_SHIFT) - Dbl.DOUBLE_EXP_OFFSET);
        val = msb & Dbl.DOUBLE_VAL_MASK;
        if (exp == Dbl.DOUBLE_EXP_MIN - 1) {
            ++exp;
        } else {
            val |= Dbl.DOUBLE_VAL_EXTRA;
        }
    }

    public void normalize() {
        if (exp < Dbl.DOUBLE_EXP_MIN) {
            toExp(Dbl.DOUBLE_EXP_MIN);
        }
        if (val == 0) {
            return;
        } 
        final int msb = msb(val);
        upBy(Math.max(msb - Dbl.DOUBLE_EXP_SHIFT, Dbl.DOUBLE_EXP_MIN - exp));
        if (val < Dbl.DOUBLE_VAL_EXTRA) {
            --exp;
        }
}

    public double toDouble() throws ValueException {
        Dbl clone = new Dbl(this);
        clone.normalize();
        if (clone.val == 0) {
            return 0;
        }
        if (clone.exp > Dbl.DOUBLE_EXP_MAX) {
            throw new ValueException(String.format("val * 2^(%d) is not finite for Dbl", this.val, this.exp));
        }
        return Double.longBitsToDouble( (neg? Dbl.DOUBLE_SIGN_MASK : 0L) | 
            ((clone.exp + Dbl.DOUBLE_EXP_OFFSET) << DOUBLE_EXP_SHIFT) | 
            (clone.val & Dbl.DOUBLE_VAL_MASK));
    }

    static final long[] BYTES = new long[7];    // for finding approx bit count quick comparison
    static {
        for (int i = 0; i < 7; ++i) {
            BYTES[i] = 1L << (8 * (i + 1));
        }
    }
    /*
     * Find the msb of (val), with lsb as 0.
     */
    static int msb(long val) {
        if (val == 0) {
            return 0;
        }
        if (val < 0) {
            return 63;
        }
        int i = 0;
        for (; i < 7; ++i) {
            if (val < BYTES[i]) {
                break;
            }
        }
        int msb = i * 8; 
        long cmp = 1L << (msb + 1);
        for (; msb < 63; ++msb, cmp <<= 1) {
            if (val < cmp) {
                break;
            }
        }
        return msb;
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
     * Round so that this.exp += shift
     *  *) When shift > 0: There is no limitation.
     *  *) When shift < 0: The val bit should be less than maxValBits (with 62 as max to allow adding two + val)
     * Return ture if exp is reached.
     */
    boolean upBy(int shift, int maxValBits) {
        if (val == 0) {
            this.exp += shift;
            return true;
        } else if (shift == 0) {
            return true;
        } else if (shift > 0) {
            this.exp += shift;
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
        } else {
            final long valExtra = 1L << (maxValBits - 1);
            if (val >= valExtra) {
                return false;
            }
            shift = -shift;
            int i = 1;
            for (; (i < 7) && (val > BYTES[i]); ++i) {
                continue;
            }
            if (i < 7) {
                if (shift < (maxValBits - 8 * i)) {
                    val <<= shift;
                    this.exp -= shift;
                    return true;
                }
                shift -= (maxValBits - 8 * i);
                val <<= (maxValBits - 8 * i);
                this.exp -= (maxValBits - 8 * i);
            }
            for (; (shift > 0) && (val < valExtra); --shift) {
                val <<= 1;
                --this.exp;
            }
            return (shift == 0);
        }
    }

    boolean upBy(int shift) {
        return upBy(shift, 62);
    }

    boolean toExp(int exp) {
        return upBy(exp - this.exp);
    }
}
