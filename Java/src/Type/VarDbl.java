package Type;

/*
 * A base class for storage type for variance arithmetic.
 * It add a taylor() method for the general Taylor expansion using the following approximation:
 *  *) M(x, 2n) = (dx)^(2n) (2n-1)!!
 *  *) M(x, 2n+1) = 0
 * It implements IReal.power() using the VarDbl.taylor().
 * 
 * This class contains the following data management private methods
 *      protected void pack(final double value, final double variance, final boolean rnd, final boolean rndr, final bound);
 *      private void pack(int exp, boolean neg, long val, boolean rnd, long var, boolean rndr, long bound);
 * The constructor assume:
 *      bias = 0
 *      linear = variance
 *      
 */
public class VarDbl implements IReal {
    public VarDbl(final double value, final double variance) throws ValueException, UncertaintyException {
        pack(value, variance, false, false, BOUND_MAX);
        bias = 0;
        linear = variance;
    }
    public VarDbl(final double value) throws ValueException {
        try {
            pack(value, Double.NaN, false, false, BOUND_MAX);
            bias = 0;
            linear = variance();
        } catch (UncertaintyException e) {
        }
    }
    public VarDbl() {
        try {
            pack(0, 0, false, false, BOUND_MAX);
            bias = 0;
            linear = 0;
        } catch (ValueException | UncertaintyException e) {
        }
    }

    public VarDbl(VarDbl other) {
        this.val = other.val;
        this.var = other.var;
        this.linear = other.linear;
        this.bias = other.bias;
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
        final Dbl dbl = new Dbl(exp(), neg(), val(), rndv());
        return dbl.toDouble();
    }

    @Override
    public double uncertainty() throws UncertaintyException {
        return Math.sqrt(variance());
    }

    public double variance() throws UncertaintyException {
        if (var() == 0) {
            return 0;
        }
        try {
            final Dbl dbl = new Dbl(exp() * 2, false, var(), rndr());
            return dbl.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException(e.getMessage());
        }
    }

    public double precSq() throws ValueException, UncertaintyException {
        final double value = value();
        return variance() / value / value;
    }

    public double ldev() {
        return Math.sqrt(linear);
    }

    public double bias() {
        return bias;
    }

    @Override
    public VarDbl negate() {
        boolean neg = !neg();
        pack(exp(), neg, val(), rndv(), var(), rndr(), bound());
        return this;
    }


    @Override
    public VarDbl shift(int bits) throws ValueException, UncertaintyException {
        if (bits == 0) {
            return this;
        }

        final Dbl dLnr = new Dbl(linear);
        dLnr.exp += 2*bits;
        final double linear;
        try {
            linear = dLnr.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException(String.format(
                "Linear variance: %s shifted by %d", IReal.format(this.linear, 2), bits
            ));
        }
        if (!Double.isFinite(linear)) {
            throw new UncertaintyException(String.format(
                "Linear variance: %s shifted by %d", IReal.format(this.linear, 2), bits
            ));
        }

        final Dbl dBias = new Dbl(bias);
        dBias.exp += bits;
        final double bias;
        try {
            bias = dBias.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException(String.format(
                "Linear variance: %s shifted by %d", IReal.format(this.bias, 2), bits
            ));
        }
        if (!Double.isFinite(bias)) {
            throw new UncertaintyException(String.format(
                "Linear variance: %s shifted by %d", IReal.format(this.bias, 2), bits
            ));
        }

        int exp = exp() + bits;
        if ((Dbl.DOUBLE_EXP_MIN <= exp) && (exp <= Dbl.DOUBLE_EXP_MAX)) {
            pack(exp, neg(), val(), rndv(), var(), rndr(), bound());
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
        pack(exp, neg(), rVal.val, rVal.rndErr, rVar.val, rVar.rndErr, bound());
        this.linear = linear;
        this.bias = bias;
        return this;
    }

    @Override
    public VarDbl add(double offset) throws ValueException, UncertaintyException {
        pack(value() + offset, variance(), rndv(), rndr(), bound());
        return this;
    }

    @Override
    public VarDbl multiply(double fold) throws ValueException, UncertaintyException {
        pack(value() * fold, variance() * fold * fold, rndv(), rndr(), bound());
        bias *= fold;
        linear *= fold * fold;
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
        final double linear = this.linear * value * value * exponent * exponent;
        throw new ValueException("NotImplementedException");
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
        final double value = value() + v.value();
        final double variance = variance() + v.variance();
        final double bias = this.bias + v.bias;
        final double linear = this.linear + v.linear;
        boolean rndErr = (neg() == v.neg())
            ? ((rndv() == v.rndv())? rndv() : false) 
            : ((rndv() != v.rndv())? rndv() : false);
        final long bound = Math.min(this.bound() >> BOUND_SHIFT, v.bound() >> BOUND_SHIFT);
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s + %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s + %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        if (!Double.isFinite(bias)) {
            throw new ValueException(String.format("%s + %s = %e: %s", 
                        IReal.format(this.bias, 3), IReal.format(v.bias, 3), bias, typeName()));
        }
        if (!Double.isFinite(linear)) {
            throw new ValueException(String.format("%s + %s = %e: %s", 
                        IReal.format(this.linear, 3), IReal.format(v.linear, 3), linear, typeName()));
        }
        pack(value, variance, rndErr, (rndr() == v.rndr())? rndr() : false, bound);
        this.bias = bias;
        this.linear = linear;
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
        final double tVal = value(), tVar = variance(), oVal = o.value(), oVar = o.variance();
        final double value = tVal * oVal;
        final double variance = tVar * oVal * oVal + oVar * tVal * tVal + tVar * oVar;
        final long bound = Math.min(this.bound(), o.bound());
        final double bias = (tVal + this.bias) * (oVal + o.bias) - value;
        final double linear = this.linear * oVal * oVal + o.linear * tVal * tVal;
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        if (!Double.isFinite(bias)) {
            throw new ValueException(String.format("%s * %s = %e: %s", 
                        IReal.format(this.bias, 3), IReal.format(o.bias, 3), bias, typeName()));
        }
        if (!Double.isFinite(linear)) {
            throw new ValueException(String.format("%s * %s = %e: %s", 
                        IReal.format(this.linear, 3), IReal.format(o.linear, 3), linear, typeName()));
        }
        pack(value, variance, 
            (rndv() == o.rndv())? rndv() : false, (rndr() == o.rndr())? rndr() : false, 
            bound);
        this.linear = linear;
        this.bias = bias;
        return this;
    }
    

    private long val;
    private long var;
    private double linear;
    private double bias;

    // encode for val
    static final int EXP_SHIFT = 53;
    static final long VAL_MASK = (1L << EXP_SHIFT) - 1;
    // encode for var
    static final int BOUND_SHIFT = 56;
    static final int VAR_RND_SHIFT = 55;
    static final int VAL_RND_SHIFT = 54;
    static final int SIGN_SHIFT = 53;
    static final int BOUND_DENOM = 32;
    static final long BOUND_MAX = (1L << (Long.SIZE - BOUND_SHIFT)) - 1;
    static final long BOUND_MASK = BOUND_MAX << BOUND_SHIFT;
    static final long VAR_RND_MASK = 1L << VAR_RND_SHIFT;
    static final long VAL_RND_MASK = 1L << VAL_RND_SHIFT;
    static final long SIGN_MASK = 1L << SIGN_SHIFT;
    static final long VAR_MASK = SIGN_MASK - 1;
 
    protected int exp()     { return (int) ((val >>> EXP_SHIFT) - Dbl.DOUBLE_EXP_OFFSET); }
    protected boolean neg() { return (var & SIGN_MASK) != 0; }
    protected long val()    { return val & VAL_MASK; }
    protected boolean rndv() { return (var & VAL_RND_MASK) != 0; }
    protected long var()    { return var & VAR_MASK; }
    protected boolean rndr() { return (var & VAR_RND_MASK) != 0; }
    protected long bound()  { return var >> BOUND_SHIFT; }

    private void pack(int exp, boolean neg, long val, boolean rnd, long var, boolean rndr, long bound) {
        this.val = ((exp + Dbl.DOUBLE_EXP_OFFSET) << EXP_SHIFT) | (val & VAL_MASK);
        this.var = (bound << BOUND_SHIFT) | (rndr? VAR_RND_MASK : 0) | (rnd? VAL_RND_MASK : 0) | (neg? SIGN_MASK : 0) | (var & VAR_MASK);
    }

    private void normalize(final Dbl dVal, final Dbl dVar, final boolean rnd, boolean rndr, long bound) throws UncertaintyException {
        if (dVar.val == 0) {
            while((dVal.val < Dbl.DOUBLE_SIGN_MASK) && (dVal.exp > Dbl.DOUBLE_EXP_MIN)) {
                dVal.val <<= 1;
                --dVal.exp;
            }
            pack(dVal.exp, dVal.neg, dVal.val, false, 0L, false, BOUND_MAX);
            return;
        }
        final Round rVal = new Round(dVal.val);
        rVal.rndErr = rnd;
        final Round rVar = new Round(dVar.val);
        rVar.rndErr = rndr();
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
        pack(dVar.exp, dVal.neg, rVal.val, rVal.rndErr, rVar.val, rVar.rndErr, bound);
    }

    protected void pack(double value, double variance, boolean rnd, boolean rndr, long bound) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%.3e~%.3e: %s()", value, Math.sqrt(variance), typeName()));
        }
        if (Double.isFinite(variance)) {
            variance = Math.abs(variance);
        } else if (Double.isNaN(variance)) {
            variance = IReal.getLSB(value);
            variance *= variance;
            if (!Double.isFinite(variance)) {
                throw new UncertaintyException(String.format("%.3e~%.3e: %s()", value, Math.sqrt(variance), typeName()));
            }
        } else {
            throw new UncertaintyException(String.format("%.3e~%.3e: %s()", value, Math.sqrt(variance), typeName()));
        }
        final Dbl dVal = new Dbl(value);
        final Dbl dVar = new Dbl(variance);
        normalize(dVal, dVar, rnd, rndr, bound);
    }
}
