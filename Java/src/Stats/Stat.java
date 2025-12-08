package Stats;

public class Stat {
	private int _count;
	private double _sum;
	private double _sum2;
	private double _min = +Double.MAX_VALUE;
	private double _max = -Double.MAX_VALUE;
	private int _minAt = 0;
	private int _maxAt = 0;

	public void clear() {
		_count = 0;
		_sum = 0;
		_sum2 = 0;
		_min = Double.MAX_VALUE;
		_max = -Double.MAX_VALUE;
	}
	
	public Stat() {
	}
	public Stat( final Stat other ) {
		this._count = other._count;
		this._sum = other._sum;
		this._sum2 = other._sum2;
		this._min = other._min;
		this._max = other._max;
	}
	
	public double min() {
		return _min;
	}
	public double max() {
		return _max;
	}
	
	public int minAt() {
		return _minAt;
	}
	public int maxAt() {
		return _maxAt;
	}
	
	public boolean accum( double value, int index ) {
		if (!Double.isFinite(value)) {
			return false;
		}
		++_count;
		_sum += value;
		_sum2 += value * value;
		if (_min > value) {
			_min = value;
            _minAt = index;
		}
		if (_max < value) {
			_max = value;
            _maxAt = index;
		}
		return true;
	}
	
	public int count() {
		return _count;
	}
	
	public double sum() {
		return _sum;
	}
	
	public double avg() {
		if (_count <= 0) {
			return Double.NaN;
		}
		return _sum / _count;
	}
	
	public double var() {
		if (_count <= 0) {
			return Double.NaN;
		}
		final double avg = avg();
		return Math.max(0, _sum2 / _count - avg * avg);
	}
	
	public double dev() {
		if (_count <= 0) {
			return Double.NaN;
		}
		return Math.sqrt( var() );
	}
}
