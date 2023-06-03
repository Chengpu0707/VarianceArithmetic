package Type;


/*
 * A tool class to round up or down
 */
class Round {
    long val;
    boolean rndErr; // True for + rounding error

    Round(long val, boolean rndErr) {
        this.val = val;
        this.rndErr = rndErr;
    }
    Round(long val) {
        this(val, false);
    }

    /*
     * shift down once, which corresponds to round up significand once.
     */
    boolean upOnce() {
        if ((val & 1) != 0) {
            if (rndErr) {
                val += 1;
                rndErr = false;
            } else {
                rndErr = true;
            }
            val >>= 1;
            return true;
        } else {
            val >>= 1;
            return false;
        }
    }

    /*
     * upBy() is eqivalent to shift * upOnce()
     * The range of negative shift is not checked, presumptively it is in the range.  If not, the result may not be expected.
     *  *) 1L << 64 == 1
     *  *) 1L << 65 == 2
     * The range of positive shift is capped at Long.SIZE
     * return true if any rounding error is generated
     */
    boolean upBy(int shift) {
        if (shift <= 0) {
            val <<= -shift;
            return false;
        }
        if (shift >= Long.SIZE) {
            final boolean ret = (val != 0);
            val = 0;
            return ret;
        }
        final long msb = 1L << shift;
        final long remain = val & (msb - 1);
        final boolean ret = (remain != 0);
        final long half = msb >> 1;
        val >>>= shift;
        if ((remain > half) || ((remain == half) && rndErr)) {
            val += 1;
            rndErr = false;
        } else if (remain > 0) {
            rndErr = true;
        }
        return ret;
    }
}


