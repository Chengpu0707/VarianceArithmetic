package Type;

/*
 * Either the value or the varance is invalid when initialize VarDbl
 */
public class InitException extends Exception {
    final double value;
    final double varinace;
    public InitException(String msg, double value, double variance) {
        super(msg);
        this.value = value;
        this.varinace = variance;
    }
}
