package Type;

/**
 * Thrown when a Taylor expansion is judged unreliable (e.g. uncertainty too
 * large relative to the value to be statistically meaningful).
 */
public class NotReliableException extends Taylor1dException {
    public NotReliableException(String msg, 
            String name, final UnionArray s1dTaylor, boolean inPrec, boolean outPrec,
            VarDbl input, VarDbl value, VarDbl variance,
            int order, VarDbl newValue, VarDbl newVariance, int monotonics) {
        super(msg, name, s1dTaylor, inPrec, outPrec,
                input, value, variance,
                order, newValue, newVariance, monotonics);
    }
}
