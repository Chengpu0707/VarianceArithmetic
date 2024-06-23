package Type;

import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;

import Type.VarDbl;

/*
 * A base class for storage type for variance arithmetic.
 * 
 * <p> This class override clone() and toString() from Object.
 * 
 * <p> This class has compareTo() but does not implement Comparable<VarDbl> to allow compareTo() throw InitException.
 * 
 * <p> Basic arithmetic operations are negate(), add(), sub(), multiply().
 */
public class VarDbl {
    private double value;
    private double variance;

    static final long DOUBLE_MAX_SIGNIFICAND = (1L << 53) - 1;
    static final long PRECISE_SIGNIFICAND_TAIL_MASK = (1L << 13) - 1;
    static final double BINDING_FOR_EQUAL = 0.67448975;
    static final double VARIANCE_THRESHOLD = 1.0 /Momentum.BINDING_FOR_TAYLOR/Momentum.BINDING_FOR_TAYLOR;
    static final double TAU = 7.18e-7;
    static final int TAYLOR_CHECK_ORDER = 20;
    static final double DEVIATION_OF_LSB = 1.0 / Math.sqrt(3);

    static public double ulp(double value) {
        return Math.ulp(value) * DEVIATION_OF_LSB;
    }
    public double ulp() {
        return ulp(value());
    }
    static public double ulp(final long value) {
        long val = Math.abs(value);
        double rounding = 0;
        boolean posi = true;
        for (; DOUBLE_MAX_SIGNIFICAND < val; val >>= 1, rounding /= 2) {
            if ((val & 1L) != 0L) {
                if (posi) {
                    rounding += 1;
                    posi = false;
                } else {
                    rounding -= 1;
                    posi = true;
                }
            }
        }
        return rounding; 
    }
    static public double ulp(final BigInteger value) {
        BigInteger val = value.abs();
        if (val.compareTo(BigInteger.valueOf(Long.MAX_VALUE)) <= 0)
            return ulp(val.longValue());
        double rounding = 0;
        boolean posi = true;
        for (int cnt = val.getLowestSetBit(); DOUBLE_MAX_SIGNIFICAND < val.doubleValue(); 
                rounding /= ((cnt < ((Long.BYTES << 3) - 1))? (1L << cnt) : Math.pow(2, cnt)), 
                val = val.shiftRight(cnt), cnt = val.getLowestSetBit()) {
            if (cnt == 0) {
                if (posi) {
                    rounding -= 1;
                    posi = false;
                }
                cnt = 1;
            }
        }
        return rounding; 
    }

    
    /*
     * Constructors
     */
    public VarDbl(final double value, double dev, boolean devAsVar) 
            throws InitException {
        final double variance = devAsVar? dev : dev * dev;
        if (!Double.isFinite(value) || !Double.isFinite(variance)) {
            throw new InitException(String.format("VarDbl value=%g~%g", value, variance), 
                                     value, variance);
        }
        this.value = value;
        this.variance = variance;
    }
    public VarDbl(final double value, double dev) 
            throws InitException {
        this(value, dev, false);
    }
    public VarDbl(final double value) 
            throws InitException {
        if (!Double.isFinite(value)) {
            throw new InitException(String.format("VarDbl value=%g", value), 
                                    value, variance);
        }
        final long val = Double.doubleToLongBits(value);
        final double dev = ((val & PRECISE_SIGNIFICAND_TAIL_MASK) == 0L)? 0 : ulp(value);
        this.value = value;
        this.variance = dev * dev;
    }
    public VarDbl(final float value) 
            throws InitException {
        if (!Double.isFinite(value)) {
            throw new InitException(String.format("VarDbl value=%g", value), 
                                    value, variance);
        }
        final long val = Float.floatToIntBits(value);
        final double dev = ((val & PRECISE_SIGNIFICAND_TAIL_MASK) == 0L)? 0 : Math.ulp(value) * DEVIATION_OF_LSB;
        this.value = value;
        this.variance = dev * dev;
    }
    public VarDbl(final long value) {
        this.value = value;
        final double dev = ulp(value);
        this.variance = dev * dev;
    }
    public VarDbl() {
        this.value = 0;
        this.variance = 0;
    }

    public VarDbl(VarDbl other) {
        this.value = other.value;
        this.variance = other.variance;
    }


    @Override
    public VarDbl clone() {
        return new VarDbl(this);
    }

	static String format( double value, int precison ) {
		if (Double.isFinite(value)) {
			if ((Math.floor(value) == Math.ceil(value)) && (Math.abs(value) < 1000)) {
				return String.format("%.0f", value); 
			} else { 
				return String.format("%." + String.format("%d", Math.abs(precison)) + "e", value);
			}
		} else if (Double.isNaN(value)) {
			return "NaN";
		} else if (Double.isInfinite(value)) {
			return (value > 0)? "+Inf" : "-Inf";
		} else {
			return "???";
		}
	}

    @Override
    public String toString() {
        final double value = value();
        final double uncertainty = uncertainty();
        if (uncertainty == 0) {
            return format(value, 3);
        } else {
            return String.format("%s~%s", format(value, 3), format(uncertainty, 1));
        }
    }

    /*
     * @ return     Comparable.compareTo()
     */
    public int compareTo(final VarDbl other, final double binding) 
            throws InitException {
        final VarDbl diff = this.minus(other);
        if (diff.variance() == 0) {
            if (diff.value() < 0)
                return -1;
            if (diff.value() > 0)
                return +1;
            return 0;
        }
        final double z = diff.value() / diff.uncertainty();
        if (Math.abs(z) <= binding)
            return 0;
        return (z < 0)? -1 : +1;
    }

    public int compareTo(final VarDbl other) 
            throws InitException {
        return compareTo(other, BINDING_FOR_EQUAL);
    }

   
    public double value() {
        return this.value;
    }

    public double uncertainty() {
        return Math.sqrt(variance());
    }

    public double variance() {
        return this.variance;
    }

    public VarDbl negate() {
        this.value = - this.value;
        return this;
    }

    /*
     * @return  inPlace? this += other : this + other;
     * 
     * @param other
     * @param inPlace
     */
    public VarDbl add(final VarDbl other, boolean inPlace) throws InitException {
        if (other == null) {
            return new VarDbl(this);
        }
        final double value = value() + other.value();
        final double variance = variance() + other.variance();
        if (!Double.isFinite(value) || !Double.isFinite(variance)) {
            throw new InitException(String.format("%s + %s = %e~%e", 
                            toString(), other.toString(), value, variance), 
                        value, variance);
        }
        if ((variance == 0) 
                && (Math.abs(value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND)
                && (Math.abs(other.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) 
                && (VarDbl.DOUBLE_MAX_SIGNIFICAND <= Math.abs(value))) {
            VarDbl sum = new VarDbl(((long )value()) + ((long) other.value()));
            if (inPlace) {
                this.value = sum.value();
                this.variance = sum.variance();
                return this;
            } else {
                return sum;
            }
        }
        if (inPlace) {
            this.value = value;
            this.variance = variance;
            return this;
        } else {
            return new VarDbl(value, variance, true);
        }
    }
    public VarDbl add(final VarDbl other) throws InitException {
        return add(other, false);
    }

    /*
     * @return  inPlace? this -= other : this - other;
     * 
     * @param other
     * @param inPlace
     */
    public VarDbl minus(final VarDbl other, boolean inPlace) throws InitException {
        return this.add(other.clone().negate(), inPlace);
    }
    public VarDbl minus(final VarDbl other) throws InitException {
        return this.minus(other, false);
    }

    /*
     * @return  inPlace? this *= other : this * other;
     * 
     * @param other
     * @param inPlace
     */
    public VarDbl multiply(final VarDbl other, boolean inPlace) throws InitException {
        final double tVal = value(), tVar = variance(), oVal = other.value(), oVar = other.variance();
        final double value = tVal * oVal;
        final double variance = tVar * oVal * oVal + oVar * tVal * tVal + tVar * oVar;
        if (!Double.isFinite(value) || !Double.isFinite(variance)) {
            throw new InitException(String.format("%s * %s = %e~%e", 
                            toString(), other.toString(), value, variance), 
                        value, variance);
        }
        if ((variance == 0) 
                && (Math.abs(value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND)
                && (Math.abs(other.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) 
                && (VarDbl.DOUBLE_MAX_SIGNIFICAND <= Math.abs(value))) {
            final BigInteger op1 = BigInteger.valueOf((long) this.value());
            final BigInteger op2 = BigInteger.valueOf((long) other.value());
            final double dev = ulp(op1.multiply(op2));
            if (inPlace) {
                this.value = value;
                this.variance = dev * dev;
                return this;
            } else {
                return new VarDbl(value, dev);
            }
        }
        if (inPlace) {
            this.value = value;
            this.variance = variance;
            return this;
        } else {
            return new VarDbl(value, variance, true);
        }
    }
    public VarDbl multiply(final VarDbl other) throws InitException {
        return multiply(other, false);
    }

    /*
     * 1d polynominal
     * 
     * @see     The paper for Variance Arithmetic on Taylor expansion convergence.
     * 
     * @return  The result of taylor expansion with this as input.
     * 
     * @param sCoeff                    Polynominal coefficients
     * @param dumpPath                  If to dump the expansion to a file
     * 
     * @exception InitException         If any item in sCoeff is not finite.
     * @exception DivergentException    If the result is not finite.
     * @exception NotReliableException  If the uncertainty of the variance is too large for its value. 
     * 
     * @raise IllegalArgumentException  If the length of sCoeff is too long, or if sCoeff contains null
     */
    VarDbl polynominal(final VarDbl[] sCoeff, final String dumpPath) 
            throws InitException, DivergentException, NotReliableException, IOException {
        if (sCoeff.length >= Momentum.MAX_FACTOR*2) {
            throw new IllegalArgumentException(String.format("coefficient length %d is too large", 
                    sCoeff.length));
        }
        final int exp = sCoeff.length - 1;
        final VarDbl[] s1dTaylor = new VarDbl[2*exp + 1];
        s1dTaylor[0] = new VarDbl(sCoeff[0]);
        for (int k = 1; k < s1dTaylor.length; ++k)
            s1dTaylor[k] = new VarDbl();
        final double[] sPow = new double[exp + 1];
        sPow[0] = 1;
        if (exp >= 1)
            sPow[1] = value();
        final int[] sTaylor = new int[2*exp + 1];
        sTaylor[0] = 1;

        for (int j = 1; j < sCoeff.length; ++j) {
            final VarDbl coeff = sCoeff[j];
            if ((coeff == null)) {
                throw new IllegalArgumentException(String.format("null coefficient at index %d/%d is too large", 
                        j, sCoeff.length));
            }
            if (!Double.isFinite(coeff.value()) || !Double.isFinite(coeff.variance())) {
                throw new InitException(String.format("polynominal coefficent at %d", j), 
                        coeff.value(), coeff.variance());
            }
            sTaylor[1] = j;
            for (int k = 2; k <= j; ++k) {
                sTaylor[k] = sTaylor[k - 1] * (j + 1 - k) / k;
                if (sPow[k] == 0)
                    sPow[k] = sPow[k - 1] * value();
            }
            for (int k = 0; k <= j; ++k) {
                s1dTaylor[k].add( coeff.multiply(new VarDbl(sTaylor[k] * sPow[j - k])), true);
            }
        }
        FileWriter fw = (dumpPath == null)? null : new FileWriter(dumpPath);
        if (fw != null) {
            fw.write("polynominal\tvalue\tuncertainty\tvariance\tBinding\tMaxOrder\n");
            fw.write(String.format("%d\t%e\t%e\t%e\t%d\t%d\nIndex\t", 
                    sCoeff.length, value(), uncertainty(), variance(),
                    Momentum.BINDING_FOR_TAYLOR, Momentum.MAX_FACTOR));
            for (int n = 0; (n < s1dTaylor.length) && (s1dTaylor[n] != null); ++n) {
                fw.write(String.format("%d\t", n));
            }
            fw.write("\nCoeff Value\t");
            for (int n = 0; n < sCoeff.length; ++n) {
                fw.write(String.format("%e\t", sCoeff[n].value()));
            }
            fw.write("\nCoeff Uncertainty\t");
            for (int n = 0; n < sCoeff.length; ++n) {
                fw.write(String.format("%e\t", sCoeff[n].uncertainty()));
            }
            fw.write("\nTaylor Value\t");
            for (int n = 0; (n < s1dTaylor.length) && (s1dTaylor[n] != null); ++n) {
                fw.write(String.format("%e\t", s1dTaylor[n].value()));
            }
            fw.write("\nTaylor Uncertainty\t");
            for (int n = 0; (n < s1dTaylor.length) && (s1dTaylor[n] != null); ++n) {
                fw.write(String.format("%e\t", s1dTaylor[n].uncertainty()));
            }
            fw.write("\n");
            fw.write("2n\tExponent Value\tExponent Variance\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty");
            fw.write("\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\n");
        }
        VarDbl value = s1dTaylor[0];
        final VarDbl variance = new VarDbl();
        final VarDbl var = new VarDbl(variance());
        final VarDbl varn = new VarDbl(var);
        int n = 2;
        for ( ; (n < Momentum.MAX_FACTOR*2) && (0 < varn.value()) && (n < s1dTaylor.length); 
                n += 2, varn.multiply(var, true)) {
            final VarDbl newValue = (s1dTaylor[n] == null)
                    ? new VarDbl() 
                    : varn.multiply(s1dTaylor[n]).multiply(new VarDbl(Momentum.factor(n)));
            final VarDbl newVariance = new VarDbl();
            for (int j = 1; j < n; ++j) {
                newVariance.add(s1dTaylor[j].multiply(s1dTaylor[n - j]).multiply(Momentum.factor(n)).multiply(varn), true);
            }
            for (int j = 2; j < n; j += 2) {
                newVariance.minus(
                    s1dTaylor[j].multiply(Momentum.factor(j)).multiply(s1dTaylor[n - j]).multiply(Momentum.factor(n - j)).multiply(varn),
                    true);
            }
            value.add(newValue, true);
            variance.add(newVariance, true);
            if (fw != null) {
                fw.write(String.format("%d\t%e\t%e\t%e\t%e\t%e\t%e", 
                         n, varn.value(), varn.variance(), value.value(), value.uncertainty(), variance.value(), variance.uncertainty()));
                fw.write(String.format("\t%e\t%e\t%e\t%e\n", 
                         newValue.value(), newValue.variance(), newVariance.value(), newVariance.uncertainty()));
            }

            if ((!Double.isFinite(value.value())) || (!Double.isFinite(value.variance()))
                    || (!Double.isFinite(variance.value())) || (!Double.isFinite(variance.variance()))) {
                if (fw != null) {
                    fw.write("DivergentException\n");
                    fw.close();
                }
                throw new DivergentException("polynominal DivergentException", 
                        s1dTaylor, false, false, this, 
                        value, variance, n, newValue, newVariance);
            }
            if (variance.variance() > variance.value() * VARIANCE_THRESHOLD) {
                if (fw != null) {
                    fw.write("NotReliableException\n");
                    fw.close();
                }
                throw new NotReliableException("polynominal NotReliableException", 
                        s1dTaylor, false, false, this,
                        value, variance, n, newValue, newVariance);
            }
        }
        if (fw != null) {
            fw.write("Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\n");
            fw.write(String.format("%e\t%e\t%e\t%e\n", 
                     value.value(), value.variance(), variance.value(), variance.variance()));
            fw.close();
        }
        final double res = variance.value() + value.variance();
        return (res == 0)? new VarDbl(value.value()) : new VarDbl(value.value(), res, true);
    }
    VarDbl polynominal(final double[] sCoeff, final String dumpPath) 
            throws InitException, DivergentException, NotReliableException, IOException {
        final VarDbl[] sVar = new VarDbl[sCoeff.length];
        for (int i = 0; i < sCoeff.length; ++i)
            sVar[i] = new VarDbl(sCoeff[i]);
        return polynominal(sVar, dumpPath);
    }
    VarDbl polynominal(final VarDbl[] sCoeff) 
            throws InitException, DivergentException, NotReliableException {
        try {
            return polynominal(sCoeff, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }
    VarDbl polynominal(final double[] sCoeff) 
            throws InitException, DivergentException, NotReliableException {
        try {
            return polynominal(sCoeff, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }


    /*
     * 1d Taylor expansion.
     * 
     * @see     The paper for Variance Arithmetic on Taylor expansion convergence.
     * 
     * @return  The result of taylor expansion with this as input.
     *          If the input is precise, the uncertainty is the floating-point uncertainty of the result.
     * 
     * @param name          The name of the Taylor expansion, for exception logging.
     * @param s1dTaylor     The Taylor expansion coefficent, with f(x) as s1dTaylor[0]. It should already contains /n!.
     * @param inPrec        If to expand by input precision
     * @param outPrec       If the variance result needs to be multiplied by s1dTaylor[0].
     * @param dumpPath                      If to dump the expansion to a file
     * @param enableStabilityTruncation     If to truncate when the expansion becomes stable.
     * 
     * @exception InitException         If any item in s1dTaylor is not finite.
     * @exception DivergentException    If the result is not finite.
     * @exception NotReliableException  If the uncertainty of the variance is too large for its value. 
     * @exceptopm NotMonotonicException If the result variance does not decrease monotonically. 
     * @exceptopm NotStableException    If after maximal order expansion, the expansion is still not stable.       
     */
    VarDbl taylor(final String name, final VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec,
            final String dumpPath, boolean enableStabilityTruncation) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException, IOException {
        if (variance() == 0)
            return new VarDbl(s1dTaylor[0]);
        FileWriter fw = (dumpPath == null)? null : new FileWriter(dumpPath);
        if (fw != null) {
            fw.write("name\tvalue\tuncertainty\tvariance\tinPrec\toutPrec\tenableStabilityTruncation\tBinding\tMaxOrder\n");
            fw.write(String.format("%s\t%e\t%e\t%e\t%b\t%b\t%b\t%d\t%d\n", 
                     name, value(), uncertainty(), variance(), inPrec, outPrec, enableStabilityTruncation,
                     Momentum.BINDING_FOR_TAYLOR, Momentum.MAX_FACTOR));
            for (int n = 0; (n < s1dTaylor.length) && (s1dTaylor[n] != null); ++n) {
                fw.write(String.format("%d\t", n));
            }
            fw.write("\n");
            for (int n = 0; (n < s1dTaylor.length) && (s1dTaylor[n] != null); ++n) {
                fw.write(String.format("%e\t", s1dTaylor[n].value()));
            }
            fw.write("\n");
            for (int n = 0; (n < s1dTaylor.length) && (s1dTaylor[n] != null); ++n) {
                fw.write(String.format("%e\t", s1dTaylor[n].uncertainty()));
            }
            fw.write("\n");
            fw.write("2n\tExponent Value\tExponent Variance\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty");
            fw.write("\tlimit\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\n");
        }
        VarDbl value = outPrec? new VarDbl(1) : new VarDbl(s1dTaylor[0]);
        final VarDbl variance = new VarDbl();
        final VarDbl var = new VarDbl(inPrec? variance() /value() /value() : variance());
        final VarDbl varn = new VarDbl(var);
        VarDbl prevVariance = null;
        int n = 2;
        for ( ; (n < Momentum.MAX_FACTOR*2) && (varn.value() > 0) && (n < s1dTaylor.length) && (s1dTaylor[n - 1] != null); 
                n += 2, varn.multiply(var, true)) {
            final VarDbl newValue = (s1dTaylor[n] == null)
                    ? new VarDbl() 
                    : varn.multiply(s1dTaylor[n]).multiply(new VarDbl(Momentum.factor(n)));
            final VarDbl newVariance = new VarDbl();
            for (int j = 1; j < n; ++j) {
                newVariance.add(s1dTaylor[j].multiply(s1dTaylor[n - j]).multiply(Momentum.factor(n)).multiply(varn), true);
            }
            for (int j = 2; j < n; j += 2) {
                newVariance.minus(
                    s1dTaylor[j].multiply(Momentum.factor(j)).multiply(s1dTaylor[n - j]).multiply(Momentum.factor(n - j)).multiply(varn),
                    true);
            }
            value.add(newValue, true);
            variance.add(newVariance, true);
            final double unc = variance.value() *TAU*TAU;
            if (fw != null) {
                fw.write(String.format("%d\t%e\t%e\t%e\t%e\t%e\t%e", 
                         n, varn.value(), varn.variance(), value.value(), value.uncertainty(), variance.value(), variance.uncertainty()));
                fw.write(String.format("\t%e\t%e\t%e\t%e\t%e\n", 
                         unc, newValue.value(), newValue.variance(), newVariance.value(), newVariance.uncertainty()));
            }

            if ((!Double.isFinite(value.value())) || (!Double.isFinite(value.variance()))
                    || (!Double.isFinite(variance.value())) || (!Double.isFinite(variance.variance()))) {
                if (fw != null) {
                    fw.write("DivergentException\n");
                    fw.close();
                }
                throw new DivergentException(name + " DivergentException", 
                        s1dTaylor, inPrec, outPrec, this, 
                        value, variance, n, newValue, newVariance);
            }
            if (variance.variance() > variance.value() * VARIANCE_THRESHOLD) {
                if (fw != null) {
                    fw.write("NotReliableException\n");
                    fw.close();
                }
                throw new NotReliableException(name + " NotReliableException", 
                        s1dTaylor, inPrec, outPrec, this,
                        value, variance, n, newValue, newVariance);
            }
            if (n >= TAYLOR_CHECK_ORDER) {
                if (Math.abs(prevVariance.value()) + unc < Math.abs(newVariance.value())) {
                    if (fw != null) {
                        fw.write("NotMonotonicException\n");
                        fw.close();
                    }
                    throw new NotMonotonicException(name + " NotMonotonicException", 
                            s1dTaylor, inPrec, outPrec, this,
                            value, variance, n, newValue, newVariance, prevVariance);
                }
                if (enableStabilityTruncation
                        && ((Math.abs(newVariance.value()) < unc) || (Math.abs(newValue.value()) < Math.ulp(value.value())))) {
                    if (fw != null) {
                        fw.write("Terminated\n");
                    }
                    break;
                }
            }
            prevVariance = newVariance;
        }
        if (enableStabilityTruncation && (Momentum.MAX_FACTOR*2 <= n)) {
            if (fw != null) {
                fw.write("NotStableException\n");
                fw.close();
            }
            throw new NotStableException(name + " NotStableException", 
                    s1dTaylor, inPrec, outPrec, this,
                    value, variance, n, prevVariance);
        }
        if (outPrec) {
            value.multiply(s1dTaylor[0], true);
            variance.multiply(s1dTaylor[0], true).multiply(s1dTaylor[0], true);
        }
        if (fw != null) {
            fw.write("Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\n");
            fw.write(String.format("%e\t%e\t%e\t%e\n", 
                     value.value(), value.variance(), variance.value(), variance.variance()));
            fw.close();
        }
       return new VarDbl(value.value(), variance.value() + value.variance(), true);
    }
    
    VarDbl taylor(final String name, VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException {
        try {
            return taylor(name, s1dTaylor, inPrec, outPrec, null, true);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return this^exponenet
     */
    public VarDbl power(double exponent, final String dumpPath) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException, IOException {
        if (exponent == 0) {
            return new VarDbl(1, 0);
        }
        if (exponent == 1) {
            return new VarDbl(this);
        }
        if ((0 < exponent) && (Math.floor(exponent) == Math.ceil(exponent))) {
            final VarDbl[] sCoeff = new VarDbl[((int) exponent) + 1];
            for (int i = 0; i < exponent; ++i) {
                sCoeff[i] = new VarDbl();
            }
            sCoeff[(int) exponent] = new VarDbl(1);
            return polynominal(sCoeff, dumpPath);
        }
        final VarDbl[] sTaylor = Taylor.power(exponent);
        sTaylor[0] = new VarDbl(Math.pow(value, exponent));
        return taylor(String.format("(%s)^%g", this, exponent), 
                sTaylor, true, false, dumpPath, true);
    }
    public VarDbl power(double exponent) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException, IOException {
        try {
            return power(exponent, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return sin(x)
     */
    public VarDbl sin(final String dumpPath) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException, IOException {
        final VarDbl[] s1dTaylor = Taylor.sin(value());
        return this.taylor(String.format("sin(%s)", this), s1dTaylor, false, false, dumpPath, true);
    }
    public VarDbl sin() 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException {
        try {
            return sin(null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return log(x)
     */
    public VarDbl log(final String dumpPath) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException, IOException {
        final VarDbl[] s1dTaylor = Taylor.log();
        return this.taylor(String.format("log(%s)", this), s1dTaylor, true, false, dumpPath, true);
    }
    public VarDbl log() 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException {
        try {
            return log(null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return exp(x)
     */
    public VarDbl exp(final String dumpPath) 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException, IOException {
        final VarDbl[] s1dTaylor = Taylor.exp();
        return this.taylor(String.format("exp(%s)", this), s1dTaylor, false, true, dumpPath, true);
    }
    public VarDbl exp() 
            throws InitException, DivergentException, NotReliableException, NotMonotonicException, NotStableException {
        try {
            return exp(null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

}
