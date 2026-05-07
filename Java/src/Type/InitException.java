package Type;

/**
 * Thrown when VarDbl construction is given invalid value or variance
 * (NaN, infinite, or negative variance).
 */
public class InitException extends Exception {
    final double value;
    final double uncertainty;
    public InitException(String msg, double value, double uncertainty) {
        super(msg);
        this.value = value;
        this.uncertainty = uncertainty;
    }
}
