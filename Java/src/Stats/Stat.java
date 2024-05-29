package Stats;

public class Stat {
	private int d_count;
	private double d_sum;
	private double d_sum2;
	private double d_min = +Double.MAX_VALUE;
	private double d_max = -Double.MAX_VALUE;

	public void clear() {
		d_count = 0;
		d_sum = 0;
		d_sum2 = 0;
		d_min = Double.MAX_VALUE;
		d_max = -Double.MAX_VALUE;
	}
	
	public Stat() {
	}
	public Stat( final Stat other ) {
		this.d_count = other.d_count;
		this.d_sum = other.d_sum;
		this.d_sum2 = other.d_sum2;
		this.d_min = other.d_min;
		this.d_max = other.d_max;
	}
	
	public double min() {
		return d_min;
	}
	public double max() {
		return d_max;
	}
	
	public boolean accum( double value ) {
		if (!Double.isFinite(value)) {
			return false;
		}
		++d_count;
		d_sum += value;
		d_sum2 += value * value;
		if (d_min > value) {
			d_min = value;
		}
		if (d_max < value) {
			d_max = value;
		}
		return true;
	}
	
	public int count() {
		return d_count;
	}
	
	public double sum() {
		return d_sum;
	}
	
	public double avg() {
		if (d_count <= 0) {
			return Double.NaN;
		}
		return d_sum / d_count;
	}
	
	public double var() {
		if (d_count <= 0) {
			return Double.NaN;
		}
		final double avg = avg();
		return Math.max(0, d_sum2 / d_count - avg * avg);
	}
	
	public double dev() {
		if (d_count <= 0) {
			return Double.NaN;
		}
		return Math.sqrt( var() );
	}
}
