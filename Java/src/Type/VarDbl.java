package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;

import Type.VarDbl;



class UnionArray {
    public UnionArray(final double[] sDbl) {
        this.sDbl = sDbl;
        this.sVar = null;
    }
    public UnionArray(final VarDbl[] sVar) {
        this.sDbl = null;
        this.sVar = sVar;
    }
    final double[] sDbl;
    final VarDbl[] sVar;
}


/*
 * A base class for storage type for variance arithmetic.
 * 
 * <p> This class override clone() and toString() from Object.
 * 
 * <p> This class has compareTo() but does not implement Comparable<VarDbl> to allow compareTo() throw InitException.
 */
public class VarDbl {
    private double value;
    private double variance;

    static final Momentum.Normal momentum = new Momentum.Normal(); // 1/sqrt(3)

    static final long DOUBLE_MAX_SIGNIFICAND = (1L << 53) - 1;
    static final long PRECISE_SIGNIFICAND_TAIL_MASK = (1L << 20) - 1;
    static final double BINDING_FOR_EQUAL = 0.67448975;
    static final double VARIANCE_THRESHOLD = 1.0 /momentum.bounding/momentum.bounding;
    static final double TAU = 7.18e-7;
    static final int MIN_MONOTONIC_COUNT = 20;
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
     * this += other
     */
    public <T> VarDbl addInPlace(final VarDbl other) throws InitException {
        if (other == null) {
            return this;
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
            VarDbl sum = new VarDbl(((long) value()) + ((long) other.value()));
            this.value = sum.value();
            this.variance = sum.variance();
            return this;
        }
        this.value = value;
        this.variance = variance;
        return this;
    }
    public VarDbl addInPlace(final double other) throws InitException {
        final double value = value() + other;
        if (!Double.isFinite(value)) {
            throw new InitException(String.format("%s + %e = %e~%e", 
                            toString(), other, value, variance), 
                        value, variance);
        }
        this.value = value;
        return this;
    }
    /*
     * this + other
     */
    public VarDbl add(final VarDbl other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.addInPlace(other);
    }
    public VarDbl add(final double other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.addInPlace(other);
    }

    /*
     * this -= other
     */
    public VarDbl minusInPlace(final VarDbl other) throws InitException {
        return this.addInPlace(other.clone().negate());
    }
    public VarDbl minusInPlace(final double other) throws InitException {
        return this.addInPlace(-other);
    }
    public VarDbl minus(final VarDbl other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.minusInPlace(other);
    }
    /*
     * this - other
     */
    public VarDbl minus(final double other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.minusInPlace(other);
    }

    /*
     * this *= other
     */
    public VarDbl multiplyInPlace(final VarDbl other) throws InitException {
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
            this.value = value;
            this.variance = dev * dev;
            return this;
        }
        this.value = value;
        this.variance = variance;
        return this;
    }
    public VarDbl multiplyInPlace(final double other) throws InitException {
        final double value = value() * other;
        final double variance = variance() * other * other;
        if (!Double.isFinite(value) || !Double.isFinite(variance)) {
            throw new InitException(String.format("%s * %e = %e~%e", 
                            toString(), other, value, variance), 
                        value, variance);
        }
        this.value = value;
        this.variance = variance;
        return this;
    }
    /*
     * this * other
     */
    public VarDbl multiply(final VarDbl other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.multiplyInPlace(other);
    }
    public VarDbl multiply(final double other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.multiplyInPlace(other);
    }

    static String INPUT_HEADER =
            "name\tvalue\tuncertainty\tvariance\tinPrec\toutPrec" +
            "\tBinding\tMaxOrder\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive";
    static String EXTENSION_HEADER = 
            "2n\tmonotonics\tExponent\tMomentum\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty" +
            "\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty";
    static String OUTPUT_HEADER =
            "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty";

    /*
    1d Taylor expansion for differential series {s1dTaylor} at {input}, with {name} for logging.
    When {inPrec} is true, calculate Taylor expnasion against the precision of {input}.
    When {outPrec} is true, the result of the Taylor expnasion is the precision.
    s1dTaylor[n] should already normalized by /n!. 
    The max order of expansion is {maxOrder}, which should not exceed momentum.Normal.MAX_ORDER

    When {checkMonotonic} is true, raise {NotMonotonicException} if 
        after full expansion, the monotonic count is still less than {MIN_MONOTONIC_COUNT}.
    It should always be True.

    When {checkStability} is true, raise {NotStableException} if 
        after full expansion, the value for the last expansion term is more than 
        TAU-fold of the expansion uncertainty.
    It should always be True.

    When {checkReliablity} is true, raise {NotReliableException} if
        the precision of the result variance is more than 1/5, in which 5 is the binding factor.
    It should always be True.

    When {checkPositive} is true, raise {NotPosive} if
        the expansion variance at any order becomes negative
    It should always be True.

    Both the result value and variance are guaranteed to be finite, otherwise
        raise {NotFiniteException} 

    Dump the expansion to {dumpPath} when it is provided.
    {dumpPath} can be read back and tested using verifyDumpFile()
     */
    private VarDbl taylor(
            final String name, final UnionArray s1dTaylor, boolean inPrec, boolean outPrec,
            final String dumpPath, final UnionArray s1dPoly, boolean checkMonotonic, boolean checkStability, boolean checkPositive,
            boolean checkReliablity) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        final int length;
        final double val;
        if (s1dTaylor.sDbl != null) {
            length = s1dTaylor.sDbl.length;
            for (int i = 0; i < length; ++i) {
                if (!Double.isFinite(s1dTaylor.sDbl[i]))
                    throw new IllegalArgumentException(String.format("Invalid Taylor coefficient [%d/%d]=%s", 
                        i, length, s1dTaylor.sDbl[i]));
            }
            val = s1dTaylor.sDbl[0];
        }
        else if (s1dTaylor.sVar != null) {
            length = s1dTaylor.sVar.length;
            for (int i = 0; i < length; ++i) {
                if (s1dTaylor.sVar == null)
                    throw new IllegalArgumentException(String.format("Null Taylor coefficient [%d/%d]", 
                        i, length));
                if (!Double.isFinite(s1dTaylor.sVar[i].value()) || !Double.isFinite(s1dTaylor.sVar[i].variance())) {
                    throw new IllegalArgumentException(String.format("Invalid Taylor coefficient [%d/%d]=%s", 
                        i, length, s1dTaylor.sDbl[i]));
                }
            }
            val = s1dTaylor.sVar[0].value();
        }
        else
            throw new IllegalArgumentException("Invalid Taylor coefficient");
        if (variance() == 0) {
            if (s1dTaylor.sDbl != null)
                return new VarDbl(val);
            else 
                return new VarDbl(s1dTaylor.sVar[0]);
        }

        FileWriter fw = (dumpPath == null)? null : new FileWriter(dumpPath);
        if (fw != null) {
            fw.write(INPUT_HEADER + "\n");
            fw.write(String.format("%s\t%e\t%e\t%e\t%b\t%b\t%.3f\t%d\t%b\t%b\t%b\t%b\n", 
                     name, value(), uncertainty(), variance(), inPrec, outPrec, 
                     momentum.bounding, momentum.maxOrder,
                     checkMonotonic, checkStability, checkReliablity, checkPositive));
            if (s1dPoly != null)
                writeList(fw, "Polynomial", s1dPoly);
            writeList(fw, "Taylor", s1dTaylor);
            fw.write(EXTENSION_HEADER + "\n");
        }

        int monotonics = 0;
        boolean monotonicPrev = true;

        VarDbl value = outPrec? new VarDbl(1) : new VarDbl((s1dTaylor.sDbl != null)? s1dTaylor.sDbl[0] : s1dTaylor.sVar[0].value());
        VarDbl variance = new VarDbl();
        final double var = inPrec? variance() /value() /value() : variance();
        double varn = var;
        VarDbl prevValue = null, prevVariance = new VarDbl();
        int n = 2;
        for ( ; (n < momentum.maxOrder) && (n < length) && Double.isFinite(varn) && (varn > 0); 
                n += 2, varn *= var) {
            final VarDbl oldValue = value, oldVariance = variance;
            final VarDbl newValue = new VarDbl(), newVariance = new VarDbl();
            String infinite = null;
            int j = 0;
            try {
                if (s1dTaylor.sDbl != null) {
                    newValue.addInPlace(s1dTaylor.sDbl[n] * varn * momentum.get(n));
                    double newVar = 0;
                    for (j = 1; j < n; ++j) {
                        newVar += s1dTaylor.sDbl[j] * s1dTaylor.sDbl[n - j] * varn *
                            (momentum.get(n) - momentum.get(j) * momentum.get(n - j));
                    }
                    newVariance.addInPlace(newVar);
                } else {
                    if (!Double.isFinite(s1dTaylor.sVar[n].value()) || !Double.isFinite(s1dTaylor.sVar[n].variance())) {
                        throw new IllegalArgumentException(String.format("Invalid Taylor coefficient [%d]=%s", n, val));
                    }
                    newValue.addInPlace(s1dTaylor.sVar[n].multiply(varn * momentum.get(n)));
                    for (j = 1; j < n; ++j) {
                        newVariance.addInPlace(s1dTaylor.sVar[j].multiply(s1dTaylor.sVar[n - j]).multiply(varn *
                            (momentum.get(n) - momentum.get(j) * momentum.get(n - j))));
                    }
                }
                value.addInPlace(newValue);
                variance.addInPlace(newVariance);
            } catch(InitException e) {
                infinite = e.getMessage();
            } catch(Throwable e) {
                if (fw != null) {
                    fw.write(String.format("Exception\t%s\t%d\t%d\t%e\n", e.getMessage(), n, j, momentum.get(n)));
                }
                throw e;
            }
            if (! Double.isFinite(value.value()))
                infinite = "value infinite";
            else if (! Double.isFinite(value.variance() + variance.value()))
                infinite = "variance infinite";
            else if (! Double.isFinite(variance.variance()))
                infinite = "variance variance infinite";
            else if (Math.abs(newVariance.value()) <= Math.abs(prevVariance.value()))
                monotonics += 1;
            else if ((monotonics > MIN_MONOTONIC_COUNT) && monotonicPrev)
                monotonicPrev = false;
            else
                monotonics = 0;
            if (infinite == null) {
                prevValue = newValue;
                prevVariance = newVariance;
            }

            if (fw != null) {
                fw.write(String.format("%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
                        n, monotonics, varn, momentum.get(n), value.value(), value.uncertainty(), variance.value(), variance.uncertainty(),
                        newValue.value(), newValue.variance(), newVariance.value(), newVariance.uncertainty()));
                fw.flush();
            }
            if (infinite != null) {
                if (checkMonotonic && (monotonics >= MIN_MONOTONIC_COUNT)) {
                    value = oldValue;
                    variance = oldVariance;
                    break;
                }
                if (fw != null) {
                    fw.write(String.format("NotFiniteException\t%d\t%s\n", n, infinite));
                    fw.close();
                }
                throw new NotFiniteException("UnionArray: " + infinite,
                        name, s1dTaylor, inPrec, outPrec,
                        this, value, variance, n, newValue, newVariance, monotonics);
            }
            if (checkPositive && (variance.value() < 0)) {
                if (fw != null) {
                    fw.write("NotPositiveException\n");
                    fw.close();
                }
                throw new NotPositiveException("negative variance",
                        name, s1dTaylor, inPrec, outPrec,
                        this, value, variance, n, newValue, newVariance, monotonics);
            }
        }

        if (checkMonotonic && (varn > 0) && (monotonics < MIN_MONOTONIC_COUNT)) {
            if (fw != null) {
                fw.write(String.format("NotMonotonicException\t%d\t%d\t%d\n", n, monotonics, MIN_MONOTONIC_COUNT));
                fw.close();
            }
            throw new NotMonotonicException(String.format("Taylor1d: MIN_MONOTONIC_COUNT=%d", MIN_MONOTONIC_COUNT),
                    name, s1dTaylor, inPrec, outPrec,
                    this, value, variance, n, prevValue, prevVariance, monotonics);
        }
        final double unc = Math.sqrt(variance.value()) *TAU;
        if (checkStability && !((Math.abs(prevValue.value()) < unc) || (Math.abs(prevValue.value()) < Math.ulp(value.value())))) {
            if (fw != null) {
                fw.write(String.format("NotStableException\t%d\t%e\t%e\t%e\n", n, unc, prevValue.value(), Math.ulp(value.value())));
                fw.close();
            }
            throw new NotStableException(String.format("Taylor1d: limit=%e", unc),
                    name, s1dTaylor, inPrec, outPrec,
                    this, value, variance, n, prevValue, prevVariance, monotonics);
        }

        if (outPrec) {
            value.multiplyInPlace(val);
            variance.multiplyInPlace(val * val);
        }
        if (! Double.isFinite(value.value())) {
            if (fw != null) {
                fw.write(String.format("NotFiniteException\t%d\t%s\t%s\n", n, "value overflow", val));
                fw.close();
            }
            throw new NotFiniteException(String.format("Taylor1d: value overflow %e", val),
                    name, s1dTaylor, inPrec, outPrec,
                    this, value, variance, n, prevValue, prevVariance, monotonics);
        }
        if (! Double.isFinite(variance.value() + value.variance())) {
            if (fw != null) {
                fw.write(String.format("NotFiniteException\t%d\t%s\t%e\n", n, "variance overflow", val));
                fw.close();
            }
            throw new NotFiniteException(String.format("Taylor1d: variance overflow %e",val),
                    name, s1dTaylor, inPrec, outPrec,
                    this, value, variance, n, prevValue, prevVariance, monotonics);
        }
        if (fw != null) {
            fw.write(OUTPUT_HEADER + "\n");
            fw.write(String.format("%e\t%e\t%e\t%e\n", 
                     value.value(), value.variance(), variance.value(), variance.variance()));
            fw.close();
        }
        return new VarDbl(value.value(), variance.value() + value.variance(), true);
    }

    private void writeList(final FileWriter fw, final String Subject, final UnionArray s1dArray) throws IOException {
        final int length;
        if (s1dArray.sVar != null) {
            length = s1dArray.sVar.length;
            fw.write("VarDbl:" + Subject);
        }
        else {
            length = s1dArray.sDbl.length;
            fw.write("Dbl:" + Subject);
        }
        for (int n = 0; n < length; ++n) {
            fw.write(String.format("\t%d", n));
        }
        fw.write("\nValue");
        if (s1dArray.sVar != null) {
            for (int n = 0; n < length; ++n) {
                fw.write(String.format("\t%e", s1dArray.sVar[n].value()));
            }
            fw.write("\nUncertainty");
            for (int n = 0; n < length; ++n) {
                fw.write(String.format("\t%e", s1dArray.sVar[n].uncertainty()));
            }
        } else {
            for (int n = 0; n < length; ++n) {
                fw.write(String.format("\t%e", s1dArray.sDbl[n]));
            }
        }
        fw.write("\n");
    }
    
    public VarDbl taylor(final String name, VarDbl[] s1dTaylor, 
                boolean inPrec, boolean outPrec, final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        return taylor(name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            dumpPath, null, true, true, true, 
            true);
    }

    public VarDbl taylor(final String name, VarDbl[] s1dTaylor, 
                boolean inPrec, boolean outPrec) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException {
        try {
            return taylor(name, s1dTaylor, inPrec, outPrec, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    public VarDbl taylor(final String name, double[] s1dTaylor, 
                boolean inPrec, boolean outPrec, final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return taylor(name, new UnionArray(s1dTaylor), inPrec, outPrec, 
            dumpPath, null, true, true, true, 
            true);
    }

    public VarDbl taylor(final String name, double[] s1dTaylor, 
                boolean inPrec, boolean outPrec) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException  {
        try {
            return taylor(name, s1dTaylor, inPrec, outPrec, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
    1d Taylor expansion for polynominal at "input" with "sCoeff".
    Allow input.value() +- input.uncertainty() to include 0
    */
    private VarDbl polynominal(final UnionArray s1dCoeff, final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        final int length;
        final UnionArray s1dTaylor;
        if (s1dCoeff.sDbl != null) {
            length = s1dCoeff.sDbl.length;
            s1dTaylor = new UnionArray(new double[length * 2 - 1]);
        }
        else if (s1dCoeff.sVar != null) {
            length = s1dCoeff.sVar.length;
            s1dTaylor = new UnionArray(new VarDbl[length * 2 - 1]);
        }
        else
            throw new IllegalArgumentException("Invalid Taylor coefficient");
        if (length * 2 >= momentum.maxOrder) {
            throw new IllegalArgumentException(String.format("coefficient length %d >= %d", 
                    length, momentum.maxOrder / 2));
        }
        final int exp = length - 1;
        if (s1dCoeff.sVar != null) {
            s1dTaylor.sVar[0] = s1dCoeff.sVar[0];
            for (int k = 1; k < s1dTaylor.sVar.length; ++k)
                s1dTaylor.sVar[k] = new VarDbl();
        } else {
            s1dTaylor.sDbl[0] = s1dCoeff.sDbl[0];
        }
        
        final double[] sPow = new double[exp + 1];
        sPow[0] = 1;
        if (exp >= 1)
            sPow[1] = value();
        final int[] sTaylor = new int[2*exp + 1];
        sTaylor[0] = 1;

        for (int j = 1; j < length; ++j) {
            VarDbl coeff = null;
            if (s1dCoeff.sVar != null) {
                coeff = s1dCoeff.sVar[j];
                if ((coeff == null)) {
                    throw new IllegalArgumentException(String.format("null coefficient at index %d/%d is too large", 
                            j, length));
                }
                if (!Double.isFinite(coeff.value()) || !Double.isFinite(coeff.variance())) {
                    throw new InitException(String.format("polynominal coefficent at %d", j), 
                            coeff.value(), coeff.variance());
                }
            }
            sTaylor[1] = j;
            for (int k = 2; k <= j; ++k) {
                sTaylor[k] = sTaylor[k - 1] * (j + 1 - k) / k;
                if (sPow[k] == 0)
                    sPow[k] = sPow[k - 1] * value();
            }
            for (int k = 0; k <= j; ++k) {
                if (coeff != null)
                    s1dTaylor.sVar[k].addInPlace(coeff.multiply(sTaylor[k] * sPow[j - k]));
                else
                    s1dTaylor.sDbl[k] += s1dCoeff.sDbl[j] * sTaylor[k] * sPow[j - k];
            }
        }
        return taylor(String.format("Poly[%d]", exp), 
                    s1dTaylor, false, false, dumpPath,
                    s1dCoeff, false, false, true, true);
    }
    VarDbl polynominal(final VarDbl[] sCoeff, final String dumpPath) 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, IOException {
        return polynominal(new UnionArray(sCoeff), dumpPath);
    }
    VarDbl polynominal(final VarDbl[] sCoeff) 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        try {
            return polynominal(sCoeff, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }
    VarDbl polynominal(final double[] sCoeff, final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        return polynominal(new UnionArray(sCoeff), dumpPath);
    }
    VarDbl polynominal(final double[] sCoeff) 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        try {
            return polynominal(sCoeff, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }


    /*
     * @return this^exponenet
     */
    public VarDbl pow(final double exponent, final String dumpPath) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        if (exponent == 0) {
            return new VarDbl(1, 0);
        }
        if (exponent == 1) {
            return new VarDbl(this);
        }
        if ((0 < exponent) && (Math.floor(exponent) == Math.ceil(exponent))) {
            final int exp = (int) exponent;
            final double[] sCoeff = new double[exp + 1];
            for (int i = 0; i < exp; ++i) {
                sCoeff[i] = 0;
            }
            sCoeff[exp] = 1;
            return polynominal(sCoeff, dumpPath);
        }
        final double[] sTaylor = new double[momentum.maxOrder];
        sTaylor[0] = 0;
        sTaylor[1] = exponent;
        double exp = exponent - 1;
        for (int i = 2; i < momentum.maxOrder; ++i, --exp) {
            sTaylor[i] = sTaylor[i - 1] * exp / i;
        }
        sTaylor[0] = Math.pow(value, exponent);
        if (!Double.isFinite(sTaylor[0]))
            throw new IllegalArgumentException(String.format("(%s)^%s=%s", this, exponent, sTaylor[0]));
        return taylor(String.format("(%s)^%s", this, exponent), 
                sTaylor, true, true, dumpPath);
    }
    public VarDbl pow(double exponent) 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        try {
            return pow(exponent, null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return sin(x)
     */
    public VarDbl sin(final String dumpPath) 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, IOException {
        final double[] sTaylor = new double[momentum.maxOrder];
        final double x = value();
        sTaylor[0] = Math.sin(x);
        double fac = 1;
        for (int i = 1; i < momentum.maxOrder; ++i, fac /= i) {
            switch (i%4) {
            case 0:
                sTaylor[i] = Math.sin(x) * fac; 
                break;
            case 1:
                sTaylor[i] = Math.cos(x) * fac;   
                break;
            case 2:
                sTaylor[i] = -Math.sin(x) * fac; 
                break;
            case 3:
                sTaylor[i] = -Math.cos(x) * fac; 
                break;
            }
        }
        return this.taylor(String.format("sin(%s)", this), sTaylor, false, false, dumpPath);
    }
    public VarDbl sin() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        try {
            return sin(null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return exp(x)
     */
    public VarDbl exp(final String dumpPath) 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, IOException {
        final double[] sTaylor = new double[momentum.maxOrder];
        sTaylor[0] = Math.exp(value());
        double fac = 1;
        for (int i = 1; i < momentum.maxOrder; ++i, fac /= i) {
            sTaylor[i] = 1.0 *fac;
        }
        return this.taylor(String.format("exp(%s)", this), sTaylor, false, true, dumpPath);
    }
    public VarDbl exp() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        try {
            return exp(null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    /*
     * @return log(x)
     */
    public VarDbl log(final String dumpPath) 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, IOException {
        final double[] sTaylor = new double[momentum.maxOrder];
        sTaylor[0] = Math.log(value());
        for (int i = 1; i < momentum.maxOrder; ++i) {
            sTaylor[i] =((i%2) == 1)? +1.0/i : -1.0/i;
        }
        return this.taylor(String.format("log(%s)", this), sTaylor, true, false, dumpPath);
    }
    public VarDbl log() 
            throws InitException, NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException {
        try {
            return log(null);
        } catch (IOException ex) {
            assert false;
            return null;
        }
    }

    static public class DumpResult {
        final Map<String, String> sInput;
        final UnionArray s1dPoly;
        final UnionArray s1dTaylor;
        final int expansions;
        final String lastLine;

        DumpResult(final Map<String, String> sInput, final UnionArray s1dPoly, final UnionArray s1dTaylor, 
                    final int expansions, final String lastLine) {
            this.sInput = sInput;
            this.s1dPoly = s1dPoly;
            this.s1dTaylor = s1dTaylor;
            this.expansions = expansions;
            this.lastLine = lastLine;
        }
    }

    static public DumpResult readDumpFile(final String dumpPath, final VarDbl res, final String excption) {
        final double REPR_DELTA = 1e-6;
        try (BufferedReader br = new BufferedReader(new FileReader(dumpPath))) {
            String line = br.readLine();
            assertEquals(line, VarDbl.INPUT_HEADER);
            String[] sHeader = line.split("\t");
            line = br.readLine();
            String[] sValue = line.split("\t");
            assertEquals(sHeader.length, sValue.length);
            final Map<String, String> sInput = new HashMap<>();
            for (int i = 0; i < sValue.length; ++i)
                sInput.put(sHeader[i], sValue[i]);

            line = br.readLine();
            sHeader = line.split("\t");
            String[] sWord = sHeader[0].split(":");
            assertTrue((sWord[0].equals("VarDbl")) || (sWord[0].equals("Dbl")));
            assertTrue((sWord[1].equals("Taylor")) || (sWord[1].equals("Polynomial")));
            final UnionArray s1dPoly;
            if (sWord[1].equals("Polynomial")) {
                s1dPoly = readList(br, sHeader);
                line = br.readLine();
                sHeader = line.split("\t");
                sWord = sHeader[0].split(":");
                assertTrue((sWord[0].equals("VarDbl")) || (sWord[0].equals("Dbl")));
                assertEquals(sWord[1], "Taylor");
            } else
                s1dPoly = null;
            final UnionArray s1dTaylor = readList(br, sHeader);

            line = br.readLine();
            assertEquals(line, VarDbl.EXTENSION_HEADER);
            sHeader = VarDbl.EXTENSION_HEADER.split("\t");
            int i = 2;
            for (; (line = br.readLine()) != null; i += 2) {
                sValue = line.split("\t");
                if (sHeader.length != sValue.length)
                    break;
                try {
                    assertEquals(Integer.parseInt(sValue[0]), i);
                } catch (NumberFormatException e) {
                    fail(String.format("Taylor expansion %d=%s parsing error %s: %s", i, sValue[0], e, line));
                }
            }
            if (excption == null) {
                assertEquals(line, VarDbl.OUTPUT_HEADER);
                sHeader = line.split("\t");
                line = br.readLine();
                sValue = line.split("\t");
                assertEquals(sHeader.length, sValue.length);
                try {
                    if (res.value() != 0)
                        assertEquals(Double.parseDouble(sValue[0]) / res.value(), 1, REPR_DELTA);
                    else
                        assertEquals(Double.parseDouble(sValue[0]), res.value(), REPR_DELTA);
                    assertEquals((Double.parseDouble(sValue[1]) + Double.parseDouble(sValue[2])) / res.variance(), 1, REPR_DELTA);
                } catch (NumberFormatException ex) {
                    fail(ex.getMessage());
                }
                } else {
                assertEquals(sValue[0], excption);
            }
            return new DumpResult(sInput, s1dPoly, s1dTaylor, i, line);
        } catch (IOException ex) {
            fail(ex.getMessage());
            return null;
        }
    }

    static private UnionArray readList(final BufferedReader br, final String[] sHeader) throws IOException {
        final String[] sWord = sHeader[0].split(":");
        assertTrue((sWord[0].equals("VarDbl")) || (sWord[0].equals("Dbl")));;
        assertTrue(2 <= sWord.length);
        for (int i = 1; i < sHeader.length; ++i) {
            try {
                assertEquals(Integer.parseInt(sHeader[i]), i - 1);
            } catch (NumberFormatException e) {
                fail(String.format("%s index %d=%s parsing error %s", sWord[1], i - 1, sHeader[i], e));
                return null;
            }
        }
        String line = br.readLine();
        final String[] sValue = line.split("\t");
        assertEquals(sHeader.length, sValue.length);
        assertEquals(sValue[0], "Value");
        final double[] sVal = new double[sHeader.length - 1];
        for (int i = 1; i < sHeader.length; ++i) {
            try {
                sVal[i - 1] = Double.parseDouble(sValue[i]);
            } catch (NumberFormatException e) {
                fail(String.format("%s value [%d]=%s parsing error %s", sWord[1], i - 1, sValue[i], e));
                return null;
            }
        }
        if (sWord[0].equals("VarDbl")) {
            line = br.readLine();
            final String[] sUnc = line.split("\t");
            assertEquals(sHeader.length, sUnc.length);
            assertEquals(sUnc[0], "Uncertainty");
            final VarDbl[] sVar = new VarDbl[sHeader.length];
            for (int i = 1; i < sHeader.length; ++i) {
                try {
                    sVar[i - 1] = new VarDbl(sVal[i -  1], Double.parseDouble(sUnc[i]));
                } catch (NumberFormatException e) {
                    fail(String.format("%s uncertiaty [%d]=%s parsing error %s", sWord[1], i - 1, sUnc[i], e));
                    return null;
                } catch (InitException e) {
                    fail(String.format("%s VarDbl [%d]=(%s,%s) parsing error %s", sWord[1], i - 1, sValue[i], sUnc[i], e));
                    return null;
                }
            }
            return new UnionArray(sVar);
        } else {
            return new UnionArray(sVal);
        }
    }

}
