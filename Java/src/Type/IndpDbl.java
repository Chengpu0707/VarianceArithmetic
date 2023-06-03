package Type;

public class IndpDbl extends VarDbl {
    public IndpDbl(double value, double variance) throws ValueException, UncertaintyException {
        super(value, variance);
    }
    public IndpDbl(long value) throws ValueException {
        super(value);
    }

    public IndpDbl() {
        super();
    }

    public IndpDbl(IndpDbl other) {
        super(other);
    }


    @Override
    public String typeName() {
         return "IndpDbl";
    }

    @Override
    public IReal power(double exponent) throws ValueException, UncertaintyException {
        if (exponent == 0) {
            return new IndpDbl(1, 0);
        }
        if (exponent == 1) {
            return new IndpDbl(this);
        }
        final double value = Math.pow(value(), exponent);
        final double variance = variance() * value * value * exponent * exponent;
        return new IndpDbl(value, variance);
    }

    @Override
    public IReal add(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (!(other instanceof IndpDbl)) {
            throw new TypeException();
        }
        final IndpDbl sum = new IndpDbl();
        return add((VarDbl) other, sum);
    }

    @Override
    public IReal multiply(final IReal other) throws TypeException, ValueException, UncertaintyException {
        if (!(other instanceof IndpDbl)) {
            throw new TypeException();
        }
        final IndpDbl prod = new IndpDbl(), ot = (IndpDbl) other;
        final double variance = ot.variance();
        prod.init(value() * other.value(), 
                  variance() * other.value() * other.value() + variance * value() * value(),
                  (rnd() == ot.rnd())? rnd() : false);
        return prod;
    }
    

    
}
