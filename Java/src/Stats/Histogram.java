package Stats;

/*
 * To construct a distogram:
 * 	*) symetric distributed between -d_maxRange and +d_maxRange
 *  *) with divids for each 1
 */
public class Histogram {

	public final int d_maxRange;
	public final int d_divids;
	protected final int d_center;
	protected final int[] d_sHistogram;
	protected final Stat d_stat = new Stat();
	protected int d_actRange = 0;
	
	public Histogram( int maxRange, int divids ) {
		d_maxRange = maxRange;
		d_divids = Math.abs(divids);
		d_center = d_divids * maxRange;
		d_sHistogram = new int[1 + d_center * 2];
	}
	public Histogram( int range ) {
		this( range, 10 );
	}
	public Histogram() {
		this( 3 );
	}

	public int divids() {
		return d_divids;
	}
	public int maxRange() {
		return d_maxRange;
	}
	
	
	public boolean accum( double value ) {
		try {
			final int idx = (int) Math.rint( value * d_divids ) + d_center;
			if ((idx < 0) || (idx >= d_sHistogram.length)) {
				return false;
			}
			d_stat.accum(value);
			++d_sHistogram[idx];
			final int actRange = Math.abs(idx - d_center);
			if (d_actRange < actRange) {
				d_actRange = actRange;
			}
			return true;
		} catch (Exception ex) {
			return false;
		}
	}
	
	public int actRange() {
		return d_actRange;
	}
	public double[] histo( int range ) {
		if ((d_stat.count() <= 0) || (d_stat.count() > Double.MAX_VALUE)) {
			return null;
		}
		if (range <= 0) {
			return new double[] {1};
		}
		final double s = 1.0/ d_stat.count();
		double[] sRes = new double[1 + 2*range];
		for (int i = 0; i < sRes.length; ++i) {
			final int idx = d_center - range + i;
			if ((idx < 0) || (idx >= d_sHistogram.length)) {
				sRes[i] = 0;
			} else {
				sRes[i] = d_sHistogram[idx] * s;
			}
		}
		return sRes;
	}
	public double[] histo() {
		return histo( d_maxRange * d_divids );
	}
	
	public Stat stat() {
		return new Stat( d_stat );
	}
}
