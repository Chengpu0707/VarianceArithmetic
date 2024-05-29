package Type;

/*
 * The Taylor expansion result has large uncetainty for either value or variance
 */
public class NotStableException extends Exception {
    public final VarDbl[] s1dTaylor;
    public final boolean inPrec, outPrec;
    public final VarDbl input, value, variance;
    public final int order;
    public final  VarDbl prevVariance;
    public NotStableException(String msg, 
                    final VarDbl[] s1dTaylor, boolean inPrec, boolean outPrec,
                    VarDbl input, VarDbl value, VarDbl variance,
                    int order, VarDbl prevVariance) {
        super(msg);
        this.s1dTaylor = s1dTaylor;
        this.inPrec = inPrec;
        this.outPrec = outPrec;
        this.input = input;
        this.value = value;
        this.variance = variance;
        this.order = order;
        this.prevVariance = prevVariance;
    }
    
}
