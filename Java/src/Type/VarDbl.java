package Type;

/*
 * A base class for storage type for variance arithmetic.
 * It add a taylor() method for the general Taylor expansion using the following approximation:
 *  *) M(x, 2n) = (dx)^(2n) (2n-1)!!
 *  *) M(x, 2n+1) = 0
 * It implements IReal.power() using the VarDbl.taylor().
 * 
 * This class contains the following data management private methods
 *      protected void init(final double value, final double variance, final boolean rnd);
 *      private void pack(int exp, boolean neg, long val, boolean rnd, long var);
 *      
 */
public class VarDbl implements IReal {
    public VarDbl(final double value, final double variance, boolean rnd) throws ValueException, UncertaintyException {
        init(value, variance, rnd);
    }
    public VarDbl(final double value, final double variance) throws ValueException, UncertaintyException {
        init(value, variance, false);
    }
    public VarDbl(final double value) throws ValueException {
        try {
            final Dbl dVal = new Dbl(value);
            final Dbl dVar = new Dbl(value);
            dVar.val = 1;
            dVar.exp *= 2;
            normalize(dVal, dVar, false);
        } catch (UncertaintyException e) {
        }
    }
    public VarDbl() {
        try {
            init(0, 0, false);
        } catch (ValueException | UncertaintyException e) {
        }
    }

    public VarDbl(VarDbl other) {
        this(other.exp(), other.neg(), other.val(), other.rnd(), other.var());
    }

    private VarDbl(int exp, boolean neg, long val, boolean rnd, long var) {
        pack(exp, neg, val, rnd, var);
    }



    @Override
    public VarDbl clone() {
        return new VarDbl(this);
    }

    @Override
    public String toString() {
        return IReal.toString(this, "~");
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof VarDbl) {
            VarDbl varDbl = (VarDbl) other;
            return (this.val == varDbl.val) && (this.var == varDbl.var);
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        return (int) this.val;
    }


    @Override
    public String typeName() {
         return "VarDbl";
    }

    @Override
    public double value() throws ValueException {
        final Dbl dbl = new Dbl(exp(), neg(), val());
        return dbl.toDouble();
    }

    public double variance() throws UncertaintyException {
        if (var() == 0) {
            return 0;
        }
        try {
            final Dbl dbl = new Dbl(exp() * 2, false, var());
            return dbl.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException(e.getMessage());
        }
    }

    @Override
    public double uncertainty() throws UncertaintyException {
        return Math.sqrt(variance());
    }

    public double precSq() throws ValueException, UncertaintyException {
        final double value = value();
        return variance() / value / value;
    }

    @Override
    public VarDbl negate() {
        boolean neg = !neg();
        pack(exp(), neg, val(), rnd(), var());
        return this;
    }

    @Override
    public VarDbl shift(int bits) throws ValueException, UncertaintyException {
        int exp = exp() + bits;
        if ((Dbl.DOUBLE_EXP_MIN <= exp) && (exp <= Dbl.DOUBLE_EXP_MAX)) {
            pack(exp, neg(), val(), rnd(), var());
            return this;
        }
        if (exp > Dbl.DOUBLE_EXP_MAX) {
            throw new ValueException(String.format("$s: shift %d overflow exp=%d val=%d var=%d", 
                        typeName(), bits, val(), var()));
        }
        final Round rVal = new Round(val());
        final Round rVar = new Round(var());
        if (exp < Dbl.DOUBLE_EXP_MIN) {
            final int shift = Dbl.DOUBLE_EXP_MIN - exp;
            if (shift >= EXP_SHIFT) {
                rVal.val = 0;
            } else {
                rVal.upBy(shift);
            }
            if (shift * 2 >= SIGN_SHIFT) {
                rVar.val = 0;
            } else {
                rVar.upBy(shift*2);
            }
            exp = Dbl.DOUBLE_EXP_MIN;
        }
        pack(exp, neg(), rVal.val, rVal.rndErr, rVar.val);
        return this;
    }

    @Override
    public VarDbl add(double offset) throws ValueException, UncertaintyException {
        init(value() + offset, variance(), rnd());
        return this;
    }

    @Override
    public VarDbl multiply(double fold) throws ValueException, UncertaintyException {
        init(value() * fold, variance() * fold * fold, rnd());
        return this;
    }

    @Override
    public VarDbl power(double exponent) throws ValueException, UncertaintyException {
        if (exponent == 0) {
            return new VarDbl(1, 0);
        }
        if (exponent == 1) {
            return new VarDbl(this);
        }
        final double value = Math.pow(value(), exponent);
        final double variance = variance() * value * value * exponent * exponent;
        return new VarDbl(value, variance);
    }

    @Override
    public VarDbl add(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            throw new TypeException(String.format("%s: %s + null", toString(), typeName()));
        }
        if (!(other instanceof VarDbl)) {
            throw new TypeException(String.format("%s: %s + %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        VarDbl v = (VarDbl) other;
        boolean rndErr = (neg() == v.neg())
            ? ((rnd() == v.rnd())? rnd() : false) 
            : ((rnd() != v.rnd())? rnd() : false);
        init(value() + v.value(), variance() + v.variance(), rndErr);
        return this;
    }

    @Override
    public VarDbl multiply(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            throw new TypeException(String.format("%s: %s * null", toString(), typeName()));
        }
        if (!(other instanceof VarDbl)) {
            throw new TypeException(String.format("%s: %s * %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        VarDbl o = (VarDbl) other;
        final double value = value(), var = variance(), oValue = o.value(), oVar = o.variance();
        init(value * oValue, 
             var * oValue * oValue + oVar * value * value + var * oVar,
             (rnd() == o.rnd())? rnd() : false);
        return this;
    }
    

    
    static final int EXP_SHIFT = 53;
    static final long VAL_MASK = (1L << EXP_SHIFT) - 1;

    static final int RND_SHIFT = 54;
    static final int SIGN_SHIFT = 53;
    static final long RND_MASK = 1L << RND_SHIFT;
    static final long SIGN_MASK = 1L << SIGN_SHIFT;
    static final long VAR_MASK = SIGN_MASK - 1;
 
    private long val;
    private long var;

    protected int exp()     { return (int) ((val >>> EXP_SHIFT) - Dbl.DOUBLE_EXP_OFFSET); }
    protected boolean neg() { return (var & SIGN_MASK) != 0; }
    protected long val()    { return val & VAL_MASK; }
    protected boolean rnd() { return (var & RND_MASK) != 0; }
    protected long var()    { return var & VAR_MASK; }

    private void pack(int exp, boolean neg, long val, boolean rnd, long var) {
        this.val = ((exp + Dbl.DOUBLE_EXP_OFFSET) << EXP_SHIFT) | (val & VAL_MASK);
        this.var = ((rnd? RND_MASK : 0) | (neg? SIGN_MASK : 0) | (var & VAR_MASK));
    }

    private void normalize(final Dbl dVal, final Dbl dVar, final boolean rnd) throws UncertaintyException {
        if (dVar.val == 0) {
            while((dVal.val < Dbl.DOUBLE_SIGN_MASK) && (dVal.exp > Dbl.DOUBLE_EXP_MIN)) {
                dVal.val <<= 1;
                --dVal.exp;
            }
            pack(dVal.exp, dVal.neg, dVal.val, false, 0);
            return;
        }
        final Round rVal = new Round(dVal.val);
        rVal.rndErr = rnd;
        final Round rVar = new Round(dVar.val);
        if ((dVar.val & Dbl.DOUBLE_VAL_EXTRA) != 0) {
            final int shift = Dbl.DOUBLE_EXP_SHIFT - SIGN_SHIFT + 1;
            rVar.upBy(shift);
            dVar.exp += shift;
        }
        if ((dVar.exp % 2) != 0) {
            rVar.upOnce();
            dVar.exp += 1;
        }
        dVar.exp /= 2;
        while (rVar.val > VAR_MASK) {
            rVar.upBy(2);
            ++dVar.exp;
        }

        int shift = dVar.exp - dVal.exp;
        if (shift > 0) {
            if (shift > Dbl.DOUBLE_EXP_SHIFT) {
                rVal.val = 0;
            } else {
                rVal.upBy(shift);
            }
        } else if (shift < 0) {
            shift = -shift;
            if (shift * 2 > SIGN_SHIFT) {
                rVar.val = 0;
            } else {
                rVar.upBy(shift * 2);
            }
            dVar.exp += shift;
        }

        // When variance > Double.MAX_EXTRA, after rouding, it become 2^53 * 2^971
        // When value = Double.MAX_EXTRA, it can never be forced to round
        dVal.val = rVar.val;
        dVal.exp = dVar.exp * 2;
        try {
            if (!Double.isFinite(dVal.toDouble())) {
                throw new UncertaintyException(String.format("Fail to normalize val=%d^%d var=%d^%d: %s", 
                            dVal.val, dVal.exp, dVar.val, dVar.exp, typeName()));
            }
        } catch (ValueException e) {
            throw new UncertaintyException(String.format("Fail to normalize val=%d^%d var=%d^%d: %s", 
                            dVal.val, dVal.exp, dVar.val, dVar.exp, typeName()));
        }
        pack(dVar.exp, dVal.neg, rVal.val, rVal.rndErr, rVar.val);
    }

    protected void init(final double value, final double variance, final boolean rnd) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("Init %.3e~%.3e: %s", value, variance, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("Init %.3e~%.3e: %s", value, variance, typeName()));
        }
        final Dbl dVal = new Dbl(value);
        final Dbl dVar = new Dbl(variance);
        normalize(dVal, dVar, rnd);
    }
}
