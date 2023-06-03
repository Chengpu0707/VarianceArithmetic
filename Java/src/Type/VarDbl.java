package Type;

import java.lang.Math;
import java.util.HashMap;
import java.util.Map;

/*
 * A 128-bit varrance representation which is consistent with a double when the variance is zero
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
            init(value, Double.NaN, false);
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
    public String typeName() {
         return "VarDbl";
    }

    @Override
    public double value() throws ValueException {
        final Dbl dbl = new Dbl(exp(), neg(), val());
        return dbl.toDouble();
    }

    @Override
    public double uncertainty() throws UncertaintyException {
        final long var = var();
        if (var == 0) {
            return 0;
        }
        try {
            final double variance = variance();
            if (variance > 0) {
                return Math.sqrt(variance);
            }
        } catch (UncertaintyException e) {
        }
        final Dbl dbl = new Dbl(exp(), false, var);
         try {
            return Math.sqrt(var) * dbl.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException();
        }
    }

    public double variance() throws UncertaintyException {
        if (var() == 0) {
            return 0;
        }
        final Dbl dbl = new Dbl(exp() * 2, false, var());
        try {
            return dbl.toDouble();
        } catch (ValueException e) {
            throw new IReal.UncertaintyException();
        }
    }

    public double bias() {
        return bias;
    }

    public double precSq() throws ValueException, UncertaintyException {
        final double value = value();
        return variance() / value / value;
    }

    @Override
    public String toString() {
        return IRealTool.toString(this);
    }

    @Override
    public IReal negate() {
        return new VarDbl(exp(), !neg(), val(), rnd(), var());
    }

    @Override
    public IReal shift(int bits) throws ValueException, UncertaintyException {
        int exp = exp() + bits;
        if ((exp <= Dbl.DOUBLE_EXP_MAX) && (exp >= Dbl.DOUBLE_EXP_MIN)) {
            return new VarDbl(exp, neg(), val(), rnd(), var());
        }
        if (exp > Dbl.DOUBLE_EXP_MAX) {
            throw new ValueException();
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
        if (rVar.val == 0) {
            if (rVal.val == 0) {
                return new VarDbl();
            }
            final Dbl dbl = new Dbl(exp, neg(), rVal.val);
            return new VarDbl(dbl.toDouble(), 0);
        }
        return new VarDbl(exp, neg(), rVal.val, rVal.rndErr, rVar.val);
    }

    @Override
    public IReal scale(double value) throws ValueException, UncertaintyException {
        VarDbl varDbl = new VarDbl();
        varDbl.init(value() * value, variance() * value * value, rnd());
        return varDbl;
    }

    protected static final int powerMaxOrder = 64;
    protected static final Map<Double, double[]> ssPowerCoeff = new HashMap<Double, double[]>();

    protected double[] powerTaylor(double exponent) {
        final double[] sTaylor = new double[ powerMaxOrder ];
        sTaylor[0] = 0;
        sTaylor[1] = exponent;
        --exponent;
        for (int i = 2; (i < powerMaxOrder) && (sTaylor[i - 1] != 0); ++i, --exponent) {
            sTaylor[i] = sTaylor[i - 1] * exponent / i;
        }
        return sTaylor;
    }

    protected double[] getPowerCoeff(double exponent) {
        if ((exponent == 0)) {
            return null;
        }
        if (ssPowerCoeff.containsKey(exponent)) {
            return ssPowerCoeff.get(exponent);
        }
        final double[] sTaylor = powerTaylor(exponent);
        final double[] sCoeff = new double[ powerMaxOrder ];
        for (int i = 1; i < powerMaxOrder; ++i) {
            for (int j = 1; j < i; ++j) {
                sCoeff[i] += sTaylor[j] * sTaylor[i - j];
            }
        }
        ssPowerCoeff.put(exponent, sCoeff);
        return sCoeff;
    }

    /*
     * Assume the input to be Gaussian at (value +/- dev)^exponent, calculate momentum of order 
     */
    protected double powerMomentum(double value, double variance, double exponent, int order) {
        if ((variance <= 0) || (exponent == 0)) {
            return Double.NaN;
        }
        if ((exponent == 1) || (order <= 0)) {
            return 0;
        }
        final double dev = Math.sqrt( variance );
        final double divid = dev / 8;
        double momentum = 0;
        double sum = 0;
        for (double x = -6*dev; x <= +6*dev; x += divid) {
            final double y = Math.pow(x, 1 / exponent);
            if (!Double.isFinite(y)) {
                continue;
            }
            final double z = (y - x) / dev;
            final double den = Math.pow(y, 1 / exponent - 1) / exponent * 
                               Math.exp(-0.5 * z * z) / dev;
            final double mmt = Math.pow(y, order) + den;
            if (!Double.isFinite(mmt)) {
                continue;
            }
            sum += den;
            momentum += mmt;
        }
        return momentum / sum;
    }

    @Override
    public IReal power( double exponent ) throws ValueException, UncertaintyException {
        if (exponent == 0) {
            return new VarDbl(1, 0);
        }
        if (exponent == 1) {
            return new VarDbl(this);
        }
        final double value = Math.pow(value(), exponent);
        if (!Double.isFinite(value)) {
            throw new ValueException();
        }
        double[] sCoeff = getPowerCoeff(exponent);
        double variance = 0;
        for (int i = 1; i < powerMaxOrder; ++i) {
            variance += sCoeff[i] * powerMomentum(value(), variance(), exponent, i);
        }
        return new VarDbl(value, variance);
    }

    @Override
    public IReal add(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (!(other instanceof VarDbl)) {
            throw new TypeException();
        }
        final VarDbl sum = new VarDbl();
        return add((VarDbl) other, sum);
    }

    protected IReal add(final VarDbl other, final VarDbl sum) throws TypeException, ValueException, UncertaintyException {
        // customized rounding to make sure that the value is accurate to the last digit
        // see TestVarDblAdd.testAddRoundingError()
        final double variance = variance() + other.variance();
        final Dbl dThis = new Dbl(exp(), neg(), val());
        final Dbl dOther = new Dbl(other.exp(), other.neg(), other.val());
        boolean rThis = rnd(), rOther = other.rnd();
        if (dThis.exp > dOther.exp) {
            final int shift = dThis.exp - dOther.exp;
            if (shift >= EXP_SHIFT) {
                dOther.val = 0;
                rOther = (other.val() != 0);
            } else {
                final Round rnd = new Round(other.val(), other.rnd());
                rnd.upBy(shift);
                dOther.val = rnd.val;
                rOther = rnd.rndErr;
            }
            dOther.exp = dThis.exp;
        } else if (dThis.exp < dOther.exp) {
            final int shift = dOther.exp - dThis.exp;
            if (shift >= EXP_SHIFT) {
                dThis.val = 0;
                rThis = (val() != 0);
            } else {
                final Round rnd = new Round(val(), rnd());
                rnd.upBy(shift);
                dThis.val = rnd.val;
                rThis = rnd.rndErr;
            }
            dThis.exp = dOther.exp;
        }
        
        if (dThis.neg == dOther.neg) {
            if (rThis == rOther) {
                sum.init(value() + other.value(), variance, rThis);
            } else {
                dThis.val += dOther.val;
                sum.init(dThis.toDouble(), variance, false);                    
            }
        } else if (rThis != rOther) {
            sum.init(value() + other.value(), variance, rThis);
        } else if (dThis.val >= dOther.val) {
            dThis.val -= dOther.val;
            sum.init(dThis.toDouble(), variance, false); 
        } else {
            dOther.val -= dThis.val;
            sum.init(dOther.toDouble(), variance, false);
        }
        return sum;
    }

    @Override
    public IReal multiply(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (!(other instanceof VarDbl)) {
            throw new TypeException();
        }
        final VarDbl prod = new VarDbl(), ot = (VarDbl) other;
        final double variance = ot.variance();
        prod.init(value() * other.value(), 
                  variance() * other.value() * other.value() + variance * value() * value() + variance() * variance,
                  (rnd() == ot.rnd())? rnd() : false);
        return prod;
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
    private double bias;

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
                throw new IReal.UncertaintyException();
            }
        } catch (ValueException e) {
            throw new IReal.UncertaintyException();
        }
        pack(dVar.exp, dVal.neg, rVal.val, rVal.rndErr, rVar.val);
    }

    protected void init(final double value, final double variance, final boolean rnd) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException();
        }
        if (Double.isFinite(variance) && Double.isNaN(variance)) {
            throw new IReal.UncertaintyException();
        }
        final Dbl dVal = new Dbl(value);
        final Dbl dVar;
        if (Double.isNaN(variance)) {
            dVar = new Dbl(value);
            dVar.val = 1;
            dVar.exp *= 2;
        } else {
            dVar = new Dbl(variance);
        }
        normalize(dVal, dVar, rnd);
    }
}
