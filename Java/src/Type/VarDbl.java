package Type;


/*
 * A base class for storage type for variance arithmetic.
 * 
 * This class contains the following data management private methods
 *      protected void pack(final double value, final double variance, final boolean rnd, final boolean rndr, final leak);
 *      private void pack(int exp, boolean neg, long val, boolean rnd, long var, boolean rndr, long leak);
 *      private void pack(final Dbl dVal, final Dbl dVar, long leak)
 *      
 */
public class VarDbl implements IReal {
    public VarDbl(final double value, final double variance) throws ValueException, UncertaintyException {
        pack(value, variance, false, false, LEAK_VAL_MAX);
    }
    public VarDbl(final double value) throws ValueException {
        try {
            pack(value, Double.NaN, false, false, LEAK_VAL_MAX);
        } catch (UncertaintyException e) {
        }
    }
    public VarDbl() {
        try {
            pack(0, 0, false, false, LEAK_VAL_MAX);
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

    public double leakage()  { 
        final long leak = (var & LEAK_MASK) >>> LEAK_SHIFT;
        long exp = leak >> LEAK_EXP_SHIFT;
        long val = leak & LEAK_VAL_MASK;
        return (val << exp) * LEAK_FACTOR; 
    }

    @Override
    public VarDbl negate() {
        boolean neg = !neg();
        pack(exp(), neg, val(), rndv(), var(), rndr(), leak());
        return this;
    }


    @Override
    public VarDbl shift(int bits) throws ValueException, UncertaintyException {
        if (bits == 0) {
            return this;
        }

        int exp = exp() + bits;
        if ((Dbl.DOUBLE_EXP_MIN <= exp) && (exp <= Dbl.DOUBLE_EXP_MAX)) {
            pack(exp, neg(), val(), rndv(), var(), rndr(), leak());
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
        pack(exp, neg(), rVal.val, rVal.rndErr, rVar.val, rVar.rndErr, leak());
        return this;
    }

    @Override
    public VarDbl add(double offset) throws ValueException, UncertaintyException {
        pack(value() + offset, variance(), rndv(), rndr(), leak());
        return this;
    }

    @Override
    public VarDbl multiply(double fold) throws ValueException, UncertaintyException {
        pack(value() * fold, variance() * fold * fold, rndv(), rndr(), leak());
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
            return this;
        }
        if (!(other instanceof VarDbl)) {
            throw new TypeException(String.format("%s: %s + %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        return add((VarDbl) other);
    }

    public VarDbl add(final VarDbl other) throws ValueException, UncertaintyException {
        if (other == null) {
            return this;
        }
        final double value = value() + other.value();
        final double variance = variance() + other.variance();
        boolean rndErr = (neg() == other.neg())
            ? ((rndv() == other.rndv())? rndv() : false) 
            : ((rndv() != other.rndv())? rndv() : false);
        final double tLkg = this.leakage(), oLkg = other.leakage();
        final long leak = leak(tLkg + oLkg - tLkg * oLkg);
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s + %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s + %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        pack(value, variance, rndErr, (rndr() == other.rndr())? rndr() : false, leak);
        return this;
    }

    @Override
    public VarDbl multiply(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            return this;
        }
        if (!(other instanceof VarDbl)) {
            throw new TypeException(String.format("%s: %s * %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        return multiply((VarDbl) other);
    }

    public VarDbl multiply(final VarDbl other) throws ValueException, UncertaintyException {
        final double tVal = value(), tVar = variance(), oVal = other.value(), oVar = other.variance();
        final double value = tVal * oVal;
        final double variance = tVar * oVal * oVal + oVar * tVal * tVal + tVar * oVar;
        final double tLkg = this.leakage(), oLkg = other.leakage();
        final long leak = leak(tLkg + oLkg - tLkg * oLkg);
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        pack(value, variance, 
            (rndv() == other.rndv())? rndv() : false, (rndr() == other.rndr())? rndr() : false, 
            leak);
        return this;
    }

    /*
     * 1d Taylor expansion.
     * 
     * @param name:         the name of the Taylor expansion, for exception logging.
     * @param s1dTaylor:    the Taylor expansion coefficent, with f(x) as s1dTaylor[0].  It should already contains /n!.
     * @param byPrec:       if to expand by precision
     * @param bounding:     the bounding factor.   
     */
    VarDbl taylor(final String name, double[] s1dTaylor, boolean byPrec, double bounding) throws ValueException, UncertaintyException {
        final int maxN = s1dTaylor.length;
        if (maxN < 2) {
            throw new ValueException(String.format("Taylor expansion with invalid coefficient of length %s", java.util.Arrays.toString(s1dTaylor)));
        }
        double value = s1dTaylor[0];
        double variance = 0;
        double var = byPrec? precSq() : variance();
        double varn = var;
        for (int n = 2; n < Momentum.maxN*2; n += 2, varn *= var) {
            value += s1dTaylor[n] * Momentum.factor(n, bounding) * varn;
            for (int j = 1; j < n; ++j) {
                variance += s1dTaylor[j] * s1dTaylor[n - j] * Momentum.factor(n, bounding) * varn;
            }
            for (int j = 2; j < n; j += 2) {
                variance -= s1dTaylor[j] * Momentum.factor(j, bounding) * 
                            s1dTaylor[n - j] * Momentum.factor(n - j, bounding) * 
                            varn;
            }
        }
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s(%s) = value %e: %s", name, toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s(%s) = variance %e: %s", name, toString(), variance, typeName()));
        }
        final double leakage = 2 - 2* Momentum.cdf(0.5 + bounding);
        final long leak = leak(leakage);
        pack(value, variance, false, false, leak);
        return this;
    }

 
    

    private long val;
    private long var;

    // encode for val
    static final int EXP_SHIFT = 53;
    static final long VAL_MASK = (1L << EXP_SHIFT) - 1;
    // encode for var
    static final int LEAK_SHIFT = 56;
    static final int VAR_RND_SHIFT = 55;
    static final int VAL_RND_SHIFT = 54;
    static final int SIGN_SHIFT = EXP_SHIFT;
    static final long VAR_MASK = VAL_MASK;
    static final long SIGN_MASK = 1L << SIGN_SHIFT;
    static final long VAR_RND_MASK = 1L << VAR_RND_SHIFT;
    static final long VAL_RND_MASK = 1L << VAL_RND_SHIFT;
    static final long LEAK_MASK = ((1L << (Long.SIZE - LEAK_SHIFT)) - 1) << LEAK_SHIFT;
    static final int LEAK_EXP_SHIFT = 4;
    static final long LEAK_VAL_MAX = (1L << LEAK_EXP_SHIFT) - 1;
    static final long LEAK_VAL_MASK = LEAK_VAL_MAX;
    static final long LEAK_EXP_MAX = (1L << (Long.SIZE - LEAK_SHIFT - LEAK_EXP_SHIFT)) - 1;
    static final long LEAK_EXP_MASK = LEAK_EXP_MAX << LEAK_EXP_SHIFT;
    static final long LEAK_EXP_OFFSET = 19;
    static final double LEAK_FACTOR = 1.0f / (1L << LEAK_EXP_OFFSET);
    static final double LEAK_MAX = (LEAK_VAL_MAX << LEAK_EXP_MAX) * LEAK_FACTOR;
 
    protected int exp()     { return (int) ((val >>> EXP_SHIFT) - Dbl.DOUBLE_EXP_OFFSET); }
    protected boolean neg() { return (var & SIGN_MASK) != 0; }
    protected long val()    { return val & VAL_MASK; }
    protected boolean rndv() { return (var & VAL_RND_MASK) != 0; }
    protected long var()    { return var & VAR_MASK; }
    protected boolean rndr() { return (var & VAR_RND_MASK) != 0; }
    protected long leak()  { 
        return var >> LEAK_SHIFT; 
    }

    static protected long leak(double leakage) throws ValueException {
        if ((leakage < 0) || (1 < leakage)) {
            throw new ValueException(String.format("Invalid leakage %g", leakage));
        } 
        if (leakage > LEAK_MAX) {
            leakage = LEAK_MAX;
        }
        Dbl dbl = new Dbl( leakage / LEAK_FACTOR );
        dbl.upBy(Dbl.DOUBLE_EXP_SHIFT - VarDbl.LEAK_EXP_SHIFT);
        if (dbl.exp() < 0) {
            dbl.upBy( -dbl.exp() );
        }
        while (VarDbl.LEAK_VAL_MAX < dbl.val()) {
            dbl.upOnce();
        }
        return (dbl.exp() << VarDbl.LEAK_EXP_SHIFT) | dbl.val();
    }

    private void pack(int exp, boolean neg, long val, boolean rnd, long var, boolean rndr, long leak) {
        this.val = ((exp + Dbl.DOUBLE_EXP_OFFSET) << EXP_SHIFT) | (val & VAL_MASK);
        this.var = (leak << LEAK_SHIFT) | (rndr? VAR_RND_MASK : 0) | (rnd? VAL_RND_MASK : 0) | (neg? SIGN_MASK : 0) | (var & VAR_MASK);
    }

    private void pack(final Dbl dVal, final Dbl dVar, long leak) throws UncertaintyException {
        if (dVar.val() == 0) {
            dVal.normalize();
            pack(dVal.exp(), dVal.neg(), dVal.val(), dVal.rndErr(), 0L, dVar.rndErr(), leak);
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

        pack(exp, dVal.neg(), dVal.val(), dVal.rndErr(), dVar.val(), dVar.rndErr(), leak);
    }

    protected void pack(double value, double variance, boolean rndv, boolean rndr, long leak) throws ValueException, UncertaintyException {
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
        pack(dVal, dVar, leak);
    }
}
