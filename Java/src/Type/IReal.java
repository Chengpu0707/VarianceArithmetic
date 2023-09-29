package Type;

/*
 * An interface to describe an imprecise value, which is a value with uncertainty.
 * 
 * Not able to specified in the interface, each derived class must have:
 *  *) A default constructor
 *  *) A copy constructor
 *  *) A constructor(double value) throws ValueException
 *  *) A constructor(double value, doulbe uncertainty) throws ValueException, UncertaintyException
 *  *) Override clone() using the copy constructor
 *  *) Override toString() using RealTool.toString(...)
 *  *) Override equals() and hashCode()
  */
public interface IReal {
	/*
	 * A static type name
	 */
	String typeName();

	/*
	 * The value and the uncertainty.  Either can be infinitive.
	 * IReal.isfinite(IReal) can be used to check the result.
	 */
	double value() throws ValueException;
	double uncertainty() throws UncertaintyException;
	
	/*
	 * Negate the value, and return this.
	 */
	IReal negate();

	/*
	 * Mutiple by 2^bits, and return this.
	 * The returned value may become 0 or infinitive, but it will not change NaN status
	 * 
	 */
	IReal shift( int bits ) throws ValueException, UncertaintyException;

	/*
	 * offset by an accurate offset, and return this.
	 * The returned value may become 0 or infinitive, but it will not change NaN status 
	 */
    IReal add( double offset ) throws ValueException, UncertaintyException;

	 /*
	 * multiple by an accurate fold, and return this.
	 * The returned value may become 0 or infinitive, but it will not change NaN status 
	 */
	IReal multiply( double scale ) throws ValueException, UncertaintyException;

	/*
	 * powered by accurate value, and return this.
	 * If exponent is false, and if 
	 */
	IReal power( double exponent ) throws ValueException, UncertaintyException;
	
	IReal add( final IReal other ) throws TypeException, ValueException, UncertaintyException;

	
	IReal multiply( final IReal other ) throws TypeException, ValueException, UncertaintyException;

	/*
	 * When the decomposed double value is no longer a double
	 */
	public static class ValueException extends Exception {
		public ValueException(String msg) {
			super(msg);
		}
	}

	/*
	 * When the variance is large than Double.MAX_VALUE
	 */
	public static class UncertaintyException extends Exception {
		public UncertaintyException(String msg) {
			super(msg);
		}
	}

	/*
	 * When two types are not compatible with each other
	 */
	public static class TypeException extends Exception {
		public TypeException(String msg) {
			super(msg);
		}
	}

	/*
	 * When the function is not implemented
	 */
	public static class NotImplementedException extends Exception {
		public NotImplementedException(String msg) {
			super(msg);
		}
	}

	/*
	 * If the range of (value, uncertainty) is finite using RealTool.isfinite(...)
	 */
	static boolean isFinite( IReal real ) {
		try {
			final double value = real.value();
			final double uncertainty = real.uncertainty();
			return Double.isFinite(value + uncertainty) && Double.isFinite(value - uncertainty);
	 	} catch (ValueException | UncertaintyException e) {
			return false;
		}
	}

	static String format( double value, int precison ) {
		if (Double.isFinite(value)) {
			if ((Math.floor(value) == Math.ceil(value)) && (Math.abs(value) < 1000)) {
				return String.format("%.0f", value); 
			} else { 
				return String.format("%." + String.format("%d", Math.abs(precison)) + "e", value);
			}
		} else if (Double.isNaN(value)) {
			return "NaN";
		} else if (Double.isInfinite(value)) {
			return (value > 0)? "+Inf" : "-Inf";
		} else {
			return "???";
		}
	}

	/*
	 * To be used to overide Object.toString()
	 */
	static String toString( IReal real, String deliminator ) {
		try {
			final double value = real.value();
			final double uncertainty = real.uncertainty();
			if (uncertainty == 0) {
				return format(value, 3);
			} else {
				return String.format("%s%s%s", format(value, 3), deliminator, format(uncertainty, 1));
			}
		} catch (ValueException e) {
			return String.format("%s:%s value", e.getMessage(), real.typeName());
		} catch (UncertaintyException e) {
			return String.format("%s:%s uncertainty", e.getMessage(), real.typeName());
		}
	}

	/*
	 * To be used when uncertainty is initiated to NaN
	 */
	static double getLSB(double value) {
		if (!Double.isFinite(value)) {
			return 0;
		}
		try {
			final Dbl d = new Dbl(value);
			final Dbl lsb = new Dbl( d.exp(), false, 1, d.rndErr());
			return lsb.toDouble();
		} catch (ValueException e) {
			return 0;
		}
	}
}

