package Type;

import Type.IReal.ValueException;

/*
 * A tool class to decompose double into different components
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
    static final long DOUBLE_VAL_NAN = 6755399441055744L;

    public Dbl( final int exp, final boolean neg, final long val ) {
      this.neg = neg;
      this.exp = exp;
      this.val = val;
    }
    public Dbl( final Dbl dbl ) {
      this(dbl.exp, dbl.neg, dbl.val);
    }

    public Dbl(final double value) {
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

    public boolean isNaN() {
      return (val == Dbl.DOUBLE_VAL_NAN) &&
             (exp == Dbl.DOUBLE_EXP_MAX + 1);
    }

    public boolean isInfinite() {
      return (val == Dbl.DOUBLE_VAL_EXTRA) &&
             (exp == Dbl.DOUBLE_EXP_MAX + 1);
    }

    public boolean isFinite() {
      return (0 <= val) && (val <= Dbl.DOUBLE_VAL_MAX) &&
             (exp <= Dbl.DOUBLE_EXP_MAX);
    }

    public double toDouble() throws ValueException {
      if (val < 0) {
        throw new IReal.ValueException();
      }
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
      while ((val < DOUBLE_VAL_EXTRA) && (exp > (Dbl.DOUBLE_EXP_MIN))) {
        val <<= 1;
        --exp;
      }
      if ((val & DOUBLE_VAL_EXTRA) == 0) {
        --exp;
      }

      if (!isFinite() && !isInfinite() && !isNaN()) {
        throw new IReal.ValueException();
      }
      return Double.longBitsToDouble( (neg? Dbl.DOUBLE_SIGN_MASK : 0L) | 
          ((exp + Dbl.DOUBLE_EXP_OFFSET) << DOUBLE_EXP_SHIFT) | 
          (val & Dbl.DOUBLE_VAL_MASK));
    }
}
