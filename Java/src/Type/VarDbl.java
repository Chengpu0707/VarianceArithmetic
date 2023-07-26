package Type;

import java.util.Arrays;

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
 *      private void pack(final Dbl dVal, final Dbl dVar, long bound)
 * 
 * The constructor assume:
 *      bias = 0
 *      linear = variance
 * The above two member are accumulated during calculation
 *      
 */
public class VarDbl implements IReal {
    public VarDbl(final double value, final double variance) throws ValueException, UncertaintyException {
        pack(value, variance, false, false, BOUND_MAX);
    }
    public VarDbl(final double value) throws ValueException {
        try {
            pack(value, Double.NaN, false, false, BOUND_MAX);
        } catch (UncertaintyException e) {
        }
    }
    public VarDbl() {
        try {
            pack(0, 0, false, false, BOUND_MAX);
        } catch (ValueException | UncertaintyException e) {
        }
    }

    public VarDbl(VarDbl other) {
        this.val = other.val;
        this.var = other.var;
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

        int exp = exp() + bits;
        if ((Dbl.DOUBLE_EXP_MIN <= exp) && (exp <= Dbl.DOUBLE_EXP_MAX)) {
            pack(exp, neg(), val(), rndv(), var(), rndr(), bound());
            return this;
        }
        if (exp > Dbl.DOUBLE_EXP_MAX) {
            throw new ValueException(String.format("%s: shift %d overflow exp=%d val=%d var=%d", 
                        typeName(), bits, exp(), val(), var()));
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
        pack(value, variance, rndErr, (rndr() == v.rndr())? rndr() : false, bound);
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
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        pack(value, variance, 
            (rndv() == o.rndv())? rndv() : false, (rndr() == o.rndr())? rndr() : false, 
            bound);
        return this;
    }

    /*
     * 1d Taylor expansion.
     * 
     * @param name:         the name of the Taylor expansion, for exception logging.
     * @param bounding:     the bounding factor.  it is infinitive if it is NaN. 
     * @param bias:         f(x + bias) - f(x)
     * @param s1dTaylor:    the Taylor expansion coefficent, with f(x) as s1dTaylor[0].  It should already contains /n!.
     * @param byPrec:       if to expand by precision
     */
    VarDbl taylor(final String name, double bounding, double bias, double[] s1dTaylor, boolean byPrec) throws ValueException, UncertaintyException {
        final int maxN = s1dTaylor.length;
        if (maxN < 2) {
            throw new ValueException(String.format("Taylor expansion with invalid coefficient of length %s", java.util.Arrays.toString(s1dTaylor)));
        }
        final double value = s1dTaylor[0];
        double variance = 0;
        for (int n = 1; n <= Momentum.maxN; ++n) {
            bias += s1dTaylor[2*n] * Momentum.factor(2*n, bounding);
            for (int j = 1; j < n; ++j) {
                variance += s1dTaylor[2*j] * s1dTaylor[2*n - 2*j] * Momentum.factor(2*n, bounding);
            }
        }
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s(%s) = value %e: %s", name, toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s(%s) = variance %e: %s", name, toString(), variance, typeName()));
        }
        final long bound = Math.min(BOUND_MAX, (long) bounding * BOUND_DENOM);
        pack(value, variance, false, false, Math.min(bound, bound()));
        return this;
    }

    

    private long val;
    private long var;

    // encode for val
    static final int EXP_SHIFT = 53;
    static final long VAL_MASK = (1L << EXP_SHIFT) - 1;
    // encode for var
    static final int BOUND_SHIFT = 56;
    static final int VAR_RND_SHIFT = 55;
    static final int VAL_RND_SHIFT = 54;
    static final int SIGN_SHIFT = EXP_SHIFT;
    static final long VAR_MASK = VAL_MASK;
    static final long SIGN_MASK = 1L << SIGN_SHIFT;
    static final long VAR_RND_MASK = 1L << VAR_RND_SHIFT;
    static final long VAL_RND_MASK = 1L << VAL_RND_SHIFT;
    static final long BOUND_MAX = (1L << (Long.SIZE - BOUND_SHIFT)) - 1;
    static final long BOUND_MASK = BOUND_MAX << BOUND_SHIFT;
    static final int BOUND_DENOM = 32;
 
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

    private void pack(final Dbl dVal, final Dbl dVar, long bound) throws UncertaintyException {
        if (dVar.val() == 0) {
            dVal.normalize();
            pack(dVal.exp(), dVal.neg(), dVal.val(), dVal.rndErr(), 0L, dVar.rndErr(), BOUND_MAX);
            return;
        }
        int shift = Math.max(0, Dbl.msb(dVar.val()) - (SIGN_SHIFT - 1));
        if ((dVar.exp() % 2) != (shift % 2)) {
            ++shift;
        }
        dVar.upBy(shift);

        shift = (dVar.exp() / 2) - dVal.exp();
        final int exp;
        if (shift > 0) {
            dVal.upBy(shift);
            exp = (dVar.exp() / 2);
        } else {
            dVar.upBy(-shift * 2);
            exp = dVal.exp();
        }

        pack(exp, dVal.neg(), dVal.val(), dVal.rndErr(), dVar.val(), dVar.rndErr(), bound);
    }

    protected void pack(double value, double variance, boolean rndv, boolean rndr, long bound) throws ValueException, UncertaintyException {
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
        final Dbl dVal = new Dbl(value, rndv);
        final Dbl dVar = new Dbl(variance, rndr);
        pack(dVal, dVar, bound);
    }
}
