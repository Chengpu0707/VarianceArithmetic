package Type;

/**
 * Thrown when a Taylor expansion produces a non-finite value or variance
 * (NaN or infinite).
 */
public class NotFiniteException extends Taylor1dException {
    public NotFiniteException(String msg, 
            String name, final UnionArray s1dTaylor, boolean inPrec, boolean outPrec,
            VarDbl input, VarDbl value, VarDbl variance,
            int order, VarDbl newValue, VarDbl newVariance, int monotonics) {
        super(msg, name, s1dTaylor, inPrec, outPrec,
                input, value, variance,
                order, newValue, newVariance, monotonics);
    }
    
}
