package Type;

/*
 * An interface to describe an imprecise value, which is a value with uncertainty
 */
public interface IReal {
	String typeName();
	double value() throws ValueException;
	double uncertainty() throws UncertaintyException;

	String toString();
	
	IReal negate();
	IReal shift( int bits ) throws ValueException, UncertaintyException;
	IReal scale( double value ) throws ValueException, UncertaintyException;
	IReal power( double exponent ) throws ValueException, UncertaintyException;
	
	IReal add( final IReal other ) throws TypeException, ValueException, UncertaintyException;
	IReal multiply( final IReal other ) throws TypeException, ValueException, UncertaintyException;

	/*
	 * When the decomposed double value is no longer a double
	 */
	public static class ValueException extends Exception {
	}
	/*
	 * When the variance is large than Double.MAX_VALUE
	 */
	public static class UncertaintyException extends Exception {
	}
	/*
	 * When two types are not compatible with each other
	 */
	public static class TypeException extends Exception {
	}
	/*
	 * When the function is not implemented
	 */
	public static class NotImplementedException extends Exception {
	}
}

class IRealTool {
	static String toString( IReal real ) {
        double value = 0, uncertainty = 0;
        try {
            value = real.value();
        } catch (IReal.ValueException e) {
            if (value == Double.POSITIVE_INFINITY) {
                return "+Inf";
            }
            if (value == Double.NEGATIVE_INFINITY) {
                return "-Inf";
            }
            if (Double.isNaN(value)) {
                return "NaN";
            }
            return "??";
        }
        try {
            uncertainty = real.uncertainty();
        } catch (IReal.UncertaintyException e) {
            if (value == Double.POSITIVE_INFINITY) {
                return String.format("%.3e~inf", value);
            }
            return String.format("%.3e~??", value);
        }
        if (uncertainty == 0) {
            return String.format("%.3e", value);
        }
        return String.format("%.3e~%.1e", value, uncertainty);
	}
}
