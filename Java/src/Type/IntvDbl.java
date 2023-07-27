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
    private double value;
    private double range;

    public IntvDbl() {
        value = 0;
        range = 0;
    }

    public IntvDbl(IntvDbl intvDbl) {
        value = intvDbl.value;
        range = intvDbl.range;
    }

    private void init( double value, double range ) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException(String.format("value %.3e for IntvDbl()", value));
        }
        if (Double.isFinite(range)) {
            range = Math.abs(range);
        } else if (Double.isNaN(range)) {
            range = IReal.getLSB(value);
        } else {
            throw new IReal.UncertaintyException(String.format("range %.3e for IntvDbl()", range));
        }      
        this.value = value;
        this.range = range;
        if (!IReal.isFinite(this)) {
            throw new IReal.UncertaintyException(String.format("%s: %s", toString(), typeName()));
        }
    }

    public IntvDbl( double value, double range ) throws ValueException, UncertaintyException {
        init(value, range);
    }

    public IntvDbl( double value ) throws ValueException {
        try {
            init(value, 0);
        } catch (UncertaintyException e) {
            // should never happens
            e.printStackTrace();
        }
    }


    @Override
    public IntvDbl clone() {
        return new IntvDbl(this);
    }

    @Override
    public String toString() {
        return IReal.toString(this, "@");
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
    public IntvDbl add(double offset) throws ValueException {
        final double value = this.value + offset;
        if (!Double.isFinite(value)) {
            throw new ValueException( String.format("%s: %s << %.3e", toString(), typeName(), offset) );
        }
        this.value = value;
        return this;
    }

    @Override
    public IntvDbl multiply(double scale) throws ValueException, UncertaintyException {
        final double value = this.value * scale;
        final double range = this.range * Math.abs(scale);
        if (!Double.isFinite(value)) {
            throw new ValueException( String.format("%s: %s * %.3e", toString(), typeName(), scale) );
        }
        if (!Double.isFinite(range)) {
            throw new UncertaintyException( String.format("%s: %s * %.3e", toString(), typeName(), scale) );
        }
        this.value = value;
        this.range = range;
        return this;
    }

    @Override
    public IntvDbl power(double exponent) throws ValueException, UncertaintyException {
        if (!Double.isFinite(exponent) || !IReal.isFinite(this)) {
            throw new ValueException(String.format("pow([%.3e, %.3e]: %s, %.3e)", 
                        value - range, value + range, typeName(), exponent));
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
            throw new UncertaintyException(String.format("pow([%.3e, %.3e]: %s, %.3e)", 
                        value - range, value + range, typeName(), exponent));
        }
        if ((exponent < 0) && ((value - range) <= 0) && (0 <= (value + range))) {
            throw new UncertaintyException(String.format("pow([%.3e, %.3e]: %s, %.3e)", 
                        value - range, value + range, typeName(), exponent));
        }
        Stat mm = new Stat();
        double d;

        d = Math.pow((value - range), exponent);
        if (!Double.isFinite(d)) {
            throw new ValueException(String.format("pow(%.3e - %.3e: %s, %.3e)", 
                        value, range, typeName(), exponent));
        }
        mm.accum(d);

        d = Math.pow((value + range), exponent);
        if (!Double.isFinite(d)) {
            throw new ValueException(String.format("pow(%.3e + %.3e: %s, %.3e)",
                        value, range, typeName(), exponent));
        }
        mm.accum(d);
 
        if ((0 < range) && ((value - range) <= 0) && (0 <= (value + range))) {
            d = Math.pow(0, exponent);
            if (!Double.isFinite(d)) {
                throw new ValueException(String.format("pow([%.3e, %.3e]: %s, %.3e)", 
                            value - range, value + range, typeName(), exponent));
            }
            mm.accum(d);
        }
        value = (mm.max() + mm.min()) * 0.5;
        range = (mm.max() - mm.min()) * 0.5;
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
        final double range = uncertainty() + other.uncertainty();
        if (!Double.isFinite(value)) {
            throw new ValueException(String.format("%s + %s: %s", toString(), other.toString(), typeName()));
        }
        if (!Double.isFinite(range)) {
            throw new UncertaintyException(String.format("%s + %s: %s", toString(), other.toString(), typeName()));
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
        final double value = other.value();
        final double range = other.uncertainty();
        Stat mm = new Stat();
        double d;
        d = (this.value - this.range)*(value - range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("(%.3e - %.3e)*(%.3e - %.3e): %s)", 
                        this.value, this.range, value, range, typeName()));
        }
        mm.accum(d);
        d = (this.value - this.range)*(value + range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("(%.3e - %.3e)*(%.3e + %.3e): %s)", 
                        this.value, this.range, value, range, typeName()));
        }
        mm.accum(d);
        d = (this.value + this.range)*(value - range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("(%.3e + %.3e)*(%.3e - %.3e): %s)", 
                        this.value, this.range, value, range, typeName()));
        }
        mm.accum(d);
        d = (this.value + this.range)*(value + range);
        if (!Double.isFinite(d)) {
            throw new UncertaintyException(String.format("(%.3e + %.3e)*(%.3e + %.3e): %s)", 
                        this.value, this.range, value, range, typeName()));
        }
        mm.accum(d);
        this.value = (mm.max() + mm.min()) * 0.5;
        this.range = (mm.max() - mm.min()) * 0.5;
        return this;
    }
}
