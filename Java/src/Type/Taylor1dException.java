package Type;

public class Taylor1dException extends Exception {
    public final String name;
    public final UnionArray s1dTaylor;
    public final boolean inPrec, outPrec;
    public final int order;
    public final VarDbl input, value, variance;
    public final VarDbl newValue, newVariance;
    public final int monotonics;
    public Taylor1dException(String msg, 
            String name, final UnionArray s1dTaylor, boolean inPrec, boolean outPrec,
            VarDbl input, VarDbl value, VarDbl variance,
            int order, VarDbl newValue, VarDbl newVariance, int monotonics) {
        super(msg);
        this.name = name;
        this.s1dTaylor = s1dTaylor;
        this.inPrec = inPrec;
        this.outPrec = outPrec;
        this.input = input;
        this.value = value;
        this.variance = variance;
        this.order = order;
        this.newValue = newValue;
        this.newVariance = newVariance;
        this.monotonics = monotonics;
    }
    
    @Override
    public String toString() {
        return String.format("%s: %s for %s at %d:" +
            "val = %s, var = %s, newVal = %s, newVar = %s, " +
            "monotonics=%d", 
            this.getClass().getName(), this.name, this.input.toString(), this.order,
            this.value.toString(), this.variance.toString(), this.newValue.toString(), this.newVariance.toString(),
            this.monotonics);
    }
}
