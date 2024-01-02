package Type;

import Stats.Stat;

/*
 * A type that implement interval arithmetic for finite double.
 * During construction or operation, IntvDbl throws:
 *  *) ValueException if value is not finite
 *  *) UncertaintyException if range is not finite. 
 *  
 * During construction, if range is NaN, the range is set on the LSB of the value without throwing UncertaintyException.
 */
public class IntvDbl implements IReal{
    static final double RANGE_OF_LSB = 1.0;
    static double getLSB(double value) {
        return Dbl.getLSB(value) * RANGE_OF_LSB;
    }
    static double getLSB(long value) {
        return Dbl.getLSB(value) * RANGE_OF_LSB;
    }
    
    private double value;
    private double range;

    public IntvDbl() {
        value = 0;
        range = 0;
    }

    public IntvDbl(final IntvDbl intvDbl) {
        value = intvDbl.value;
        range = intvDbl.range;
    }

    public IntvDbl( double value, double range ) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) 
            throw new IReal.ValueException(String.format("value %.3e for IntvDbl()", value));
        if (!Double.isFinite(range)) 
            throw new IReal.UncertaintyException(String.format("range %.3e for IntvDbl()", range));
        range = Math.abs(range);
        this.value = value;
        this.range = range;
        if (!isFinite()) {
            throw new IReal.UncertaintyException(String.format("%s: %s", toString(), typeName()));
        }
    }

    public IntvDbl( double value ) throws ValueException, UncertaintyException {
        this( value, getLSB(value));
    }

    public IntvDbl(final long value) {
        // (long) this.double is casted to Long.MAX_VALUE, so the following does not work
        // if ((long) this.value == value)
        this.value = value;
        this.range = getLSB(value);
    }

    @Override
    public IntvDbl clone() {
        return new IntvDbl(this);
    }

    @Override
    public String toString() {
        return IReal.toString(this, "+-");
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof IntvDbl) {
            IntvDbl intv = (IntvDbl) other;
            return (Double.compare(this.value, intv.value) == 0) &&
                   (Double.compare(this.range, intv.range) == 0);
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        return Double.hashCode(value);
    }


    @Override
    public String typeName() {
        return "IntvDbl";
    }

    @Override
    public double value() {
        return value;
    }

    @Override
    public double uncertainty() {
        return range;
    }

    @Override
    public IntvDbl negate() {
        value = -value;
        return this;
    }

    @Override
    public IntvDbl shift(int bits) throws ValueException, UncertaintyException {
        Dbl val = new Dbl(value);
        final Dbl unc = new Dbl(range);
        val.shift(bits);
        unc.shift(bits);
        final double value, range;
        try {
            value = val.toDouble();
        } catch (ValueException e) {
            throw new ValueException(String.format("%s: %s << %d", toString(), typeName(), bits));
        }
        try {
            range = unc.toDouble();
        } catch (ValueException e) {
            throw new UncertaintyException(String.format("%s: %s << %d", toString(), typeName(), bits));
        }
        this.value = value;
        this.range = range;
        return this;
    }

    @Override
    public IntvDbl power(double exponent) throws ValueException, UncertaintyException {
        if (!Double.isFinite(exponent) || !isFinite()) {
            throw new ValueException(String.format("%s^%.3e", 
                        toString(), exponent));
        }
        if (exponent == 0) {
            value = 1;
            range = 0;
            return this;
        }
        if (exponent == 1) {
            return this;
        }
        if ((Math.floor(exponent) != Math.ceil(exponent)) && ((value - range) < 0)) {
            throw new UncertaintyException(String.format("%s^%.3e", 
                        toString(), exponent));
        }
        if ((exponent < 0) && ((value - range) <= 0) && (0 <= (value + range))) {
            throw new UncertaintyException(String.format("%s^%.3e", 
                        toString(), exponent));
        }
        Stat mm = new Stat();
        double d;

        d = Math.pow((value - range), exponent);
        if (!Double.isFinite(d)) {
            throw new ValueException(String.format("%s^%.3e: %.3e^%.3e=%.3e", 
                        toString(), exponent, value - range, exponent, d));
        }
        mm.accum(d);

        d = Math.pow((value + range), exponent);
        if (!Double.isFinite(d)) {
            throw new ValueException(String.format("%s^%.3e: %.3e^%.3e=%.3e", 
                        toString(), exponent, value + range, exponent, d));
        }
        mm.accum(d);
 
        if (Math.floor(exponent) != Math.ceil(exponent)) {
            d = Math.pow((value - range), (exponent - Dbl.getLSB(exponent)));
            if (!Double.isFinite(d)) {
                throw new ValueException(String.format("%s^%.3e: %.3e^%.3e=%.3e", 
                            toString(), exponent, value - range, exponent - Dbl.getLSB(exponent), d));
            }
            mm.accum(d);

            d = Math.pow((value + range), (exponent - Dbl.getLSB(exponent)));
            if (!Double.isFinite(d)) {
                throw new ValueException(String.format("%s^%.3e: %.3e^%.3e=%.3e", 
                            toString(), exponent, value + range, exponent - Dbl.getLSB(exponent), d));
            }
            mm.accum(d);

            d = Math.pow((value - range), (exponent + Dbl.getLSB(exponent)));
            if (!Double.isFinite(d)) {
                throw new ValueException(String.format("%s^%.3e: %.3e^%.3e=%.3e", 
                            toString(), exponent, value - range, exponent + Dbl.getLSB(exponent), d));
            }
            mm.accum(d);

            d = Math.pow((value + range), (exponent + Dbl.getLSB(exponent)));
            if (!Double.isFinite(d)) {
                throw new ValueException(String.format("%s^%.3e: %.3e^%.3e=%.3e", 
                            toString(), exponent, value + range, exponent + Dbl.getLSB(exponent), d));
            }
            mm.accum(d);
        }

        if ((0 < range) && ((value - range) <= 0) && (0 <= (value + range))) {
            d = Math.pow(0, exponent);
            if (!Double.isFinite(d)) {
                throw new ValueException(String.format("pow([%.3e, %.3e]: %s, %.3e)", 
                            value - range, value + range, typeName(), exponent));
            }
            mm.accum(d);
        }

        value = (mm.max() + mm.min()) * 0.5;
        range = (mm.max() - mm.min()) * 0.5 + Dbl.getLSB(value);
        return this;
    }

    @Override
    public IntvDbl add(IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            return this;
        }
        if (!(other instanceof IntvDbl)) {
            throw new TypeException(String.format("%s: %s + %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        return add( (IntvDbl) other);
    }

    public IntvDbl add(IntvDbl other) throws ValueException, UncertaintyException {
        if (other == null) {
            return this;
        }
        final double value = value() + other.value();
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s + %s: %s", toString(), other.toString(), typeName()));
        }
        double range = uncertainty() + other.uncertainty();
        if (!Double.isFinite(range)) {
            throw new UncertaintyException(String.format("%s + %s: %s", toString(), other.toString(), typeName()));
        }
        if (range == 0) {
            range = getLSB((long) value() + (long) other.value());
        } 
        this.value = value;
        this.range = range;
        return this;
    }

    @Override
    public IntvDbl multiply(IReal other) throws TypeException, ValueException, UncertaintyException {
        if (other == null) {
            return this;
        }
        if (!(other instanceof IntvDbl)) {
            throw new TypeException(String.format("%s: %s * %s: %s", 
                        toString(), typeName(), other.toString(), other.typeName()));
        }
        return multiply((IntvDbl) other);
    }

    public IntvDbl multiply(IntvDbl other) throws ValueException, UncertaintyException {
        if (other == null) {
            return this;
        }
        if ((this.range == 0) && (other.range == 0)) {
            final double value = this.value * other.value();
            if (!Double.isFinite(value)) 
                throw new UncertaintyException(String.format("%.3e*%.3e:%.3e %s)", 
                            this.value, other.value, value, typeName()));
            this.value = value;
            this.range = getLSB(value);
            return this;
        }
        final double value = other.value;
        final double range = other.range;
        Stat mm = new Stat();
        double d;
        d = (this.value - this.range)*(value - range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("%s * %s: (%.3e-%.3e)*(%.3e-%.3e)=%.3e)", 
                        toString(), other.toString(), this.value, this.range, value, range, d));
        }
        mm.accum(d);
        d = (this.value - this.range)*(value + range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("%s * %s: (%.3e-%.3e)*(%.3e+%.3e)=%.3e)", 
                        toString(), other.toString(), this.value, this.range, value, range, d));
        }
        mm.accum(d);
        d = (this.value + this.range)*(value - range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("%s * %s: (%.3e+%.3e)*(%.3e-%.3e)=%.3e)", 
                        toString(), other.toString(), this.value, this.range, value, range, d));
        }
        mm.accum(d);
        d = (this.value + this.range)*(value + range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("%s * %s: (%.3e+%.3e)*(%.3e+%.3e)=%.3e)", 
                        toString(), other.toString(), this.value, this.range, value, range, d));
        }
        mm.accum(d);
        this.value = mm.max() * 0.5 + mm.min() * 0.5;
        this.range = mm.max() * 0.5 - mm.min() * 0.5 + Dbl.getLSB(this.value);
        return this;
    }
}
