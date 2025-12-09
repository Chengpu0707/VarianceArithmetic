package Stats;

/*
 * To construct a distogram:
 * 	*) symetric distributed between -_maxRange and +_maxRange
 *  *) with divids for each 1
 */
public class Histogram {

	public final double _maxRange;
	public final int _divids;
	protected final int _center;
	protected final int[] _sHistogram;
	protected final Stat _stat = new Stat();
	protected int _lower = 0;
    protected int _upper = 0;
	
	public Histogram( double maxRange, int divids ) {
		_maxRange = maxRange;
		_divids = Math.abs(divids);
		_center = (int) (_divids * maxRange);
		_sHistogram = new int[1 + _center * 2];
	}
	public Histogram( int range ) {
		this( range, 5 );
	}
	public Histogram() {
		this( 3 );
	}

	public int divids() {
		return _divids;
	}
	public double maxRange() {
		return _maxRange;
	}
	
	
	public boolean accum( double value, int index ) {
		try {
			_stat.accum(value, index);
			final int idx = (int) Math.rint( value * _divids ) + _center;
			if (idx < 0) {
                ++_lower;
				return false;
			} else if (idx >= _sHistogram.length) {
                ++_upper;
                return false;
            }
			++_sHistogram[idx];
			return true;
		} catch (Exception ex) {
			return false;
		}
	}
	
	public int lower() {
		return _lower;
	}
	public int upper() {
		return _upper;
	}


	public double[] histo() {
        final int range = (int) (_maxRange * _divids);
        final int cnt = _stat.count() - lower() - upper();
		if (cnt <= 0) {
			return null;
		}
		final double s = 1.0/ cnt;
		double[] sRes = new double[1 + 2*range];
		for (int i = 0; i < sRes.length; ++i) {
			final int idx = _center - range + i;
			if ((idx < 0) || (idx >= _sHistogram.length)) {
				sRes[i] = 0;
			} else {
				sRes[i] = _sHistogram[idx] * s;
			}
		}
		return sRes;
	}
	
	public Stat stat() {
		return new Stat( _stat );
	}
}
