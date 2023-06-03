package Type;

public class IntvDbl implements IReal{

    public IntvDbl( double value, double range ) throws ValueException, UncertaintyException {
        init( value, range );
    }
    public IntvDbl( double value ) throws ValueException, UncertaintyException {
        init( value, Double.NaN );
    }
    public IntvDbl() {
        try {
            init( 0, 0 );
        } catch (ValueException | UncertaintyException e) {
        }
    }

    public IntvDbl(IntvDbl intvDbl) {
        try {
            init( intvDbl.value, intvDbl.range );
        } catch (ValueException | UncertaintyException e) {
        }
    }

    @Override
    public String typeName() {
        return "IntvDbl";
    }

    @Override
    public double value() throws ValueException {
        return value;
    }

    @Override
    public double uncertainty() throws UncertaintyException {
        return range;
    }

    @Override
    public String toString() {
        return IRealTool.toString(this);
    }

    @Override
    public IReal negate() {
        try {
            return new IntvDbl( -value, range );
        } catch (ValueException | UncertaintyException e) {
            return null;
        }
    }

    @Override
    public IReal shift(int bits) throws ValueException, UncertaintyException {
        final Dbl val = new Dbl(value());
        final Dbl unc = new Dbl(uncertainty());
        val.exp += bits;
        unc.exp += bits;
        final double value = val.toDouble(), uncertainty;
        try {
            uncertainty = unc.toDouble();
        } catch (ValueException e) {
            throw new IReal.UncertaintyException();
        }
        return new IntvDbl( value, uncertainty );
    }

    @Override
    public IReal scale(double fold) throws ValueException, UncertaintyException {
        final double value = value() * fold;
        final double range = uncertainty() * fold;
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException();
        }
        if (!Double.isFinite(range)) {
            throw new IReal.UncertaintyException();
        }
        return new IntvDbl(value, range);
    }

    @Override
    public IReal power(double exponent) throws ValueException, UncertaintyException {
        if (exponent == 0) {
            return new IntvDbl(1, 0);
        }
        if (exponent == 1) {
            return new IntvDbl(this);
        }
        if ((exponent < 0) && ((value - range) <= 0) && ((value + range) >= 0)) {
            throw new UncertaintyException();
        }
        if ((value - range) >= 0) {
            final double min = (exponent > 0)? Math.pow(value - range, exponent) : Math.pow(value + range, exponent);
            final double max = (exponent > 0)? Math.pow(value + range, exponent) : Math.pow(value - range, exponent);
            return new IntvDbl( (max + min) * 0.5, (max - min) * 0.5);
        }
        final double min_, max_;
        try {
            min_ = Math.pow((value - range), exponent);
            max_ = Math.pow((value + range), exponent);
        } catch (Throwable e) {
            throw new ValueException();
        }
        if (Double.isNaN(min_) || Double.isNaN(max_)) {
            throw new ValueException();
        }
        final double max = Math.max(min_, max_);
        double min = Math.min(min_, max_);
        if ((min > 0) && ((value - range) <= 0) && ((value + range) >= 0)) {
            min = 0;
        }
        return new IntvDbl( (max + min) * 0.5, (max - min) * 0.5);
    }

    @Override
    public IReal add(IReal other) throws TypeException, ValueException, UncertaintyException {
        if (!(other instanceof IntvDbl)) {
            throw new TypeException();
        }
        final double value = value() + other.value();
        final double range = uncertainty() + other.uncertainty();
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException();
        }
        if (!Double.isFinite(range)) {
            throw new IReal.UncertaintyException();
        }
        return new IntvDbl(value, range);
    }

    @Override
    public IReal multiply(IReal other) throws TypeException, ValueException, UncertaintyException {
        final double value = other.value();
        final double range = other.uncertainty();
        final double[] sLimt = new double[] {
            (this.value - this.range)*(value - range),
            (this.value + this.range)*(value - range),
            (this.value + this.range)*(value - range),
            (this.value + this.range)*(value + range),
        };
        double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
        for (double v : sLimt) {
            if (min > v) {
                min = v;
            }
            if (max < v) {
                max = v;
            }
        }
        if (!Double.isFinite(max) && !Double.isFinite(min)) {
            throw new IReal.ValueException();
        }
        if (!Double.isFinite(max) || !Double.isFinite(min)) {
            throw new IReal.UncertaintyException();
        }
        return new IntvDbl((max + min)/2.0, (max - min)/2.0);
    }

    private void init( double value, double range ) throws ValueException, UncertaintyException {
        if (!Double.isFinite(value)) {
            throw new IReal.ValueException();
        }
        if (Double.isFinite(range)) {
            range = Math.abs(range);
        } else if (Double.isNaN(range)) {
            final Dbl d = new Dbl(value);
            d.val = 1;
            d.neg = false;
            range = d.toDouble();
        } else {
            throw new IReal.UncertaintyException();
        }
        
        if (!Double.isFinite(value + range)){
            throw new IReal.UncertaintyException();
        }
        if (!Double.isFinite(value - range)){
            throw new IReal.UncertaintyException();
        }
        this.value = value;
        this.range = range;
    }


    private double value;
    private double range;
    
}
