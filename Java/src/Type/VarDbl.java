package Type;

import java.math.BigInteger;

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
public class VarDbl implements Comparable<VarDbl> {
    private double value;
    private double uncertainty;

    public double value() {
        return this.value;
    }

    public double uncertainty() {
        return this.uncertainty;
    }

    public double variance() {
        return this.uncertainty * this.uncertainty;
    }

    static final long DOUBLE_MAX_SIGNIFICAND = (1L << 53) - 1;
    static final long PRECISE_SIGNIFICAND_TAIL_MASK = (1L << 20) - 1;
    static final double BINDING_FOR_EQUAL = 0.67448975;
    static final double LEAKAGE_5 = 5.733031437360481e-07;
    static final double DEVIATION_OF_LSB = 1.0 / Math.sqrt(3);

    static public double ulp(double value) {
        return Math.ulp(value) * DEVIATION_OF_LSB;
    }
    static public double ulp(float value) {
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
    public VarDbl(final double value, double uncertainty) 
            throws InitException {
        if (!Double.isFinite(value) || !Double.isFinite(uncertainty)) {
            throw new InitException(String.format("VarDbl value=%g~%g", value, uncertainty), 
                                     value, uncertainty);
        }
        this.value = value;
        this.uncertainty = uncertainty;
    }
    public VarDbl(final double value) 
            throws InitException {
        if (!Double.isFinite(value)) {
            throw new InitException(String.format("VarDbl value=%g", value), 
                                    value, uncertainty);
        }
        final long val = Double.doubleToLongBits(value);
        final double uncertainty = ((val & PRECISE_SIGNIFICAND_TAIL_MASK) == 0L)? 0 : ulp(value);
        this.value = value;
        this.uncertainty = uncertainty;
    }
    public VarDbl(final float value) 
            throws InitException {
        if (!Double.isFinite(value)) {
            throw new InitException(String.format("VarDbl value=%g", value), 
                                    value, 0);
        }
        final long val = Float.floatToIntBits(value);
        final double uncertainty = ((val & PRECISE_SIGNIFICAND_TAIL_MASK) == 0L)? 0 : Math.ulp(value) * DEVIATION_OF_LSB;
        this.value = value;
        this.uncertainty = uncertainty;
    }
    public VarDbl(final long value) {
        this.value = value;
        this.uncertainty = ulp(value);
    }
    public VarDbl() {
        this.value = 0;
        this.uncertainty = 0;
    }

    public VarDbl(VarDbl other) {
        this.value = other.value;
        this.uncertainty = other.uncertainty;
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

    @Override
    public int compareTo(final VarDbl other) {
        final double val = this.value() - other.value();
        final double unc = Math.sqrt(this.variance() + other.variance());
        if (!Double.isFinite(unc))
            return 0;
         if (unc == 0)
            return (val == 0)? 0 : (val < 0)? -1 : +1;
        final double z = val / unc;
        if (Math.abs(z) <= BINDING_FOR_EQUAL)
            return 0;
        return (z < 0)? -1 : +1;
    }
    public int compareTo(final double other) throws InitException {
        return compareTo(new VarDbl(other));
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
        final double value = this.value() + other.value();
        final double uncertainty = (this.uncertainty() == 0)? other.uncertainty() :
                (other.uncertainty() == 0)? this.uncertainty() : Math.sqrt(this.variance() + other.variance());
        if (!Double.isFinite(value) || !Double.isFinite(uncertainty)) {
            throw new InitException(String.format("%s + %s = %e~%e", 
                            toString(), other.toString(), value, uncertainty), 
                        value, uncertainty);
        }
        this.value = value;
        this.uncertainty = uncertainty;
        return this;
    }
    public VarDbl addInPlace(final double other) throws InitException {
        return addInPlace(new VarDbl(other));
    }
    public VarDbl add(final VarDbl other) throws InitException {
        return new VarDbl(this).addInPlace(other);
    }
    public VarDbl add(final double other) throws InitException {
        return new VarDbl(this).addInPlace(other);
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
         return new VarDbl(this).minusInPlace(other);
    }
    public VarDbl minus(final double other) throws InitException {
        final VarDbl var = new VarDbl(this);
        return var.minusInPlace(other);
    }

    /*
     * this *= other
     */
    public VarDbl multiplyInPlace(final VarDbl other) throws InitException {
        final double value = this.value() * other.value();
        double uncertainty;
        if ((other.uncertainty() == 0) && (this.uncertainty() == 0)) {  
            final BigInteger op1 = BigInteger.valueOf((long) this.value());
            final BigInteger op2 = BigInteger.valueOf((long) other.value());
            uncertainty = ulp(op1.multiply(op2));
        } else if (other.uncertainty() == 0) {
            uncertainty = this.uncertainty() *  Math.abs(other.value());
        } else if (this.uncertainty() == 0) {
            uncertainty = other.uncertainty() * Math.abs(this.value());
        } else {
            uncertainty = Math.sqrt(
                    this.variance() * other.value() * other.value() +
                    other.variance() * this.value() * this.value() +
                    this.variance() * other.variance());
            if (uncertainty == 0) {
                uncertainty = 
                    this.uncertainty() * Math.abs(other.value()) +
                    other.uncertainty() *  Math.abs(this.value()) +
                    this.uncertainty() * other.uncertainty();
            }
        }
        if (!Double.isFinite(value) || !Double.isFinite(uncertainty)) {
            throw new InitException(String.format("%s * %s = %e~%e", 
                            toString(), other.toString(), value, uncertainty), 
                        value, uncertainty);
        }
        this.value = value;
        this.uncertainty = uncertainty;
        return this;
    }
    public VarDbl multiplyInPlace(final double other) throws InitException {
        return this.multiplyInPlace(new VarDbl(other));
    }
    public VarDbl multiply(final VarDbl other) throws InitException {
        return new VarDbl(this).multiplyInPlace(other);
    }
    public VarDbl multiply(final double other) throws InitException {
        return this.multiply(new VarDbl(other));
    }

}
