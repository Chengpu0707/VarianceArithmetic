package Type;

/*
 * The Taylor expansion result has large uncetainty for either value or variance
 */
public class NotStableException extends Taylor1dException {
    public NotStableException(String msg, 
            String name, final UnionArray s1dTaylor, boolean inPrec, boolean outPrec,
            VarDbl input, VarDbl value, VarDbl variance,
            int order, VarDbl newValue, VarDbl newVariance, int monotonics) {
        super(msg, name, s1dTaylor, inPrec, outPrec,
                input, value, variance,
                order, newValue, newVariance, monotonics);
    }
}
