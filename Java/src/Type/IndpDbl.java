package Type;

public class IndpDbl implements IReal {
    private double value;
    private double variance;

    public IndpDbl() {
        value = 0;
        variance = 0;
    }

    public IndpDbl(IndpDbl other) {
        value = other.value;
        variance = other.variance;
    }

    private void init( double value, double variance ) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException(String.format("value %.3e for IndpDbl()", value));
        }
        if (Double.isFinite(variance)) {
            variance = Math.abs(variance);
        } else if (Double.isNaN(variance)) {
            variance = IReal.getLSB(value);
            variance *= variance;
        } else {
            throw new IReal.UncertaintyException(String.format("variance %.3e for IndpDbl()", variance));
        }      
        this.value = value;
        this.variance = variance;
        if (!IReal.isFinite(this)) {
            throw new IReal.UncertaintyException(String.format("%s: %s", toString(), typeName()));
        }
    }

    public IndpDbl( double value, double variance ) throws ValueException, UncertaintyException {
        init(value, variance);
    }

    public IndpDbl( double value ) throws ValueException {
        try {
            init(value, 0);
        } catch (UncertaintyException e) {
            // should never happens
            e.printStackTrace();
        }
    }


    @Override
    public String typeName() {
         return "IndpDbl";
    }

    @Override
    public IndpDbl clone() {
        return new IndpDbl(this);
    }

    @Override
    public String toString() {
        return IReal.toString(this, "&");
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof IndpDbl) {
            IndpDbl intv = (IndpDbl) other;
            return (Double.compare(this.value, intv.value) == 0) &&
                   (Double.compare(this.variance, intv.variance) == 0);
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        return Double.hashCode(value);
    }

    @Override
    public double value() {
        return value;
    }

    @Override
    public double uncertainty() {
        return Math.sqrt(variance);
    }

    @Override
    public IndpDbl negate() {
        value = -value;
        return this;
    }

    @Override
    public IndpDbl shift(int bits) throws ValueException, UncertaintyException {
        Dbl val = new Dbl(value);
        final Dbl unc = new Dbl(variance);
        val.exp += bits;
        unc.exp += bits;
        final double value, variance;
        try {
            value = val.toDouble();
        } catch (ValueException e) {
            throw new ValueException(String.format("%s: %s << %d", toString(), typeName(), bits));
        }
        try {
            variance = unc.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException(String.format("%s: %s << %d", toString(), typeName(), bits));
        }
        this.value = value;
        this.variance = variance;
        return this;
    }

    @Override
    public IndpDbl add(double offset) throws ValueException {
        final double value = this.value + offset;
        if (!Double.isFinite(value)) {
            throw new ValueException( String.format("%s: %s << %.3e", toString(), typeName(), offset) );
        }
        this.value = value;
        return this;
    }

    @Override
    public IndpDbl multiply(double scale) throws ValueException, UncertaintyException {
        final double value = this.value * scale;
        final double variance = this.variance * Math.abs(scale);
        if (!Double.isFinite(value)) {
            throw new ValueException( String.format("%s: %s * %.3e", toString(), typeName(), scale) );
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException( String.format("%s: %s * %.3e", toString(), typeName(), scale) );
        }
        this.value = value;
        this.variance = variance;
        return this;
    }

    @Override
    public IndpDbl power(double exponent) throws ValueException, UncertaintyException {
        if (exponent == 0) {
            value = 1;
            variance = 0;
            return this;
        }
        if (exponent == 1) {
            return this;
        }
        final double value = Math.pow(this.value, exponent);
        final double f = value / this.value * exponent;
        final double variance = this.variance * f * f;
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s^%s: %s", toString(), IReal.format(exponent, 3), typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s^%s: %s", toString(), IReal.format(exponent, 3), typeName()));
        }
        this.value = value;
        this.variance = variance;
        return this;
    }

    @Override
    public IReal add(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            throw new TypeException(String.format("%s: %s + null", toString(), typeName()));
        }
        if (!(other instanceof IndpDbl)) {
            throw new TypeException(String.format("%s: %s + %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        final double value = value() + other.value();
        final double variance = uncertainty() + other.uncertainty();
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s + %s: %s", toString(), other.toString(), typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s + %s: %s", toString(), other.toString(), typeName()));
        }
        this.value = value;
        this.variance = variance;
        return this;
    }

    @Override
    public IndpDbl multiply(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            throw new TypeException(String.format("%s: %s * null", toString(), typeName()));
        }
        if (!(other instanceof IndpDbl)) {
            throw new TypeException(String.format("%s: %s * %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        final double value = this.value() * other.value();
        final double variance = this.uncertainty() * other.value() * other.value() + 
                                other.uncertainty() * this.value() * this.value();
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s*%s: %s)", toString(), other.toString(), typeName()));
        }
        if (!Double.isFinite(variance)) {
            throw new UncertaintyException(String.format("%s*%s: %s)", toString(), other.toString(), typeName()));
        }
        this.value = value;
        this.variance = variance;
        return this;
    }
   
}
