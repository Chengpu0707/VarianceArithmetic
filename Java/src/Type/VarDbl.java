package Type;


/*
 * A base class for storage type for variance arithmetic      
 */
public class VarDbl implements IReal {
    static final double BINDING_FOR_EQUAL = 0.67448975;

    static final double DEVIATION_OF_LSB = 1.0 / Math.sqrt(3);
    static public double getLSB(double value) {
        return Math.ulp(value) * DEVIATION_OF_LSB;
    }
    static public double getLSB(long value) {
        return Dbl.getLSB(value) * DEVIATION_OF_LSB;
    }

    private double value;
    private double variance;

    public VarDbl(final double value, double dev) throws ValueException, UncertaintyException {
        final double variance = dev * dev;
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s value=%f~%f", typeName(), value, variance));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s variance=%f~%f", typeName(), value, variance));
        }
        this.value = value;
        this.variance = variance;
    }
    public VarDbl(final double value) throws ValueException {
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s value=%f", typeName(), value));
        }
        this.value = value;
        final double dev = getLSB(value);
        this.variance = dev * dev;
    }
    public VarDbl(final long value) {
        this.value = value;
        final double dev = getLSB(value);
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

    @Override
    public String toString() {
        return IReal.toString(this, "~");
    }

    public int compareTo(final VarDbl other, final double binding) 
            throws ValueException, UncertaintyException {
        final VarDbl diff = clone();
        diff.minus(other);
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
            throws ValueException, UncertaintyException {
        return compareTo(other, BINDING_FOR_EQUAL);
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof VarDbl) {
            return (this.value() == ((VarDbl) other).value()) && (this.variance() == ((VarDbl) other).variance());
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        return (int) value();
    }


    @Override
    public String typeName() {
         return "VarDbl";
    }

    @Override
    public double value() {
        return this.value;
    }

    @Override
    public double uncertainty() {
        return Math.sqrt(variance());
    }

    public double variance() {
        return this.variance;
    }

    public double precSq() {
        return variance() / value() / value();
    }

    @Override
    public VarDbl negate() {
        this.value = - this.value;
        return this;
    }


    @Override
    public VarDbl shift(int bits) throws ValueException, UncertaintyException {
        if (bits == 0) {
            return this;
        }
        Dbl dbl = new Dbl(this.value);
        dbl.shift(bits);
        this.value = dbl.toDouble();

        dbl = new Dbl(this.variance);
        dbl.shift(bits * 2);
        try {
            this.variance = dbl.toDouble();
        } catch (ValueException e) {
            throw new IReal.UncertaintyException(e.getMessage());
        }
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
        final double[] sTaylor = Taylor.power(exponent);
        sTaylor[0] = Math.pow(value, exponent);
        return taylor("power", sTaylor, true, false);
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
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s + %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        double variance = variance() + other.variance();
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s + %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        if (variance == 0) {
            final double dev = getLSB((long) value() + (long) other.value());
            variance = dev * dev;
        } 
        this.value = value;
        this.variance = variance;
        return this;
    }

    public VarDbl minus(final VarDbl other) throws ValueException, UncertaintyException {
        return this.add(other.clone().negate());
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
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), value, typeName()));
        }
        double variance = tVar * oVal * oVal + oVar * tVal * tVal + tVar * oVar;
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s * %s = %e: %s", 
                        toString(), other.toString(), variance, typeName()));
        }
        if (variance == 0) {
            final int msb = Dbl.msb(Math.abs((long) value())) + Dbl.msb(Math.abs((long) other.value()));
            if (msb >= Dbl.DOUBLE_EXP_SHIFT) {
                final double dev = getLSB(value);
                variance = dev * dev;
            }
        }
        this.value = value;
        this.variance = variance;
        return this;
    }

    /*
     * 1d Taylor expansion.
     * 
     * @param name:         the name of the Taylor expansion, for exception logging.
     * @param s1dTaylor:    the Taylor expansion coefficent, with f(x) as s1dTaylor[0].  
     *                      It should already contains /n!.
     * @param inPrec:       if to expand by input precision
     * @param outPrec:      if the variance result needs to be multiplied by s1dTaylor[0]
     * @param bounding:     the bounding factor. 
     * 
     * @return:  
     */
    VarDbl taylor(final String name, double[] s1dTaylor, boolean inPrec, boolean outPrec) 
            throws ValueException, UncertaintyException {
        if (variance() == 0) {
            throw new UncertaintyException(String.format("%s(%s) = 0 variance: %s", name, toString(), typeName()));
        }
        double value = outPrec? 1 : s1dTaylor[0];
        double variance = 0;
        double var = inPrec? precSq() : variance();
        double varn = var;
        for (int n = 2; n < Momentum.maxN*2; n += 2, varn *= var) {
            value += s1dTaylor[n] * Momentum.factor(n) * varn;
            for (int j = 1; j < n; ++j) {
                variance += s1dTaylor[j] * s1dTaylor[n - j] * Momentum.factor(n) * varn;
            }
            for (int j = 2; j < n; j += 2) {
                variance -= s1dTaylor[j] * Momentum.factor(j) * 
                            s1dTaylor[n - j] * Momentum.factor(n - j) * 
                            varn;
            }
        }
        if (outPrec) {
            value *= s1dTaylor[0];
            variance *= s1dTaylor[0] * s1dTaylor[0];
        }
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s(%s) = value %e: %s", name, toString(), value, typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s(%s) = variance %e: %s", name, toString(), variance, typeName()));
        }
        this.value = value;
        this.variance = variance;
        return this;
    }
    
}
