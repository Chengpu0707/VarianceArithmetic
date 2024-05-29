package Type;

/*
 * The Taylor expansion result in infinitive value or variance
 */
public class DivergentException extends Exception {
    public final VarDbl[] s1dTaylor;
    public final boolean inPrec, outPrec;
    public final VarDbl input, value, variance;
    public final int order;
    public final  VarDbl newValue, newVariance;
    public DivergentException(String msg, 
                    final VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec,
                    VarDbl input, VarDbl value, VarDbl variance,
                    int order, VarDbl newValue, VarDbl newVariance) {
        super(msg);
        this.s1dTaylor = s1dTaylor;
        this.inPrec = inPrec;
        this.outPrec = outPrec;
        this.input = input;
        this.value = value;
        this.variance = variance;
        this.order = order;
        this.newValue = newValue;
        this.newVariance = newVariance;
    }
    
}
