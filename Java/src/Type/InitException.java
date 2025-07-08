package Type;

/*
 * Either the value or the varance is invalid when initialize VarDbl
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
