package Type;


/**
 * Abstract base for bounded moment tables used by Taylor expansion: stores the
 * truncated moment array, leakage and maxOrder for a given bounding factor.
 * Subclassed by {@link NormalMoment} for Gaussian distributions.
 *
 * @since 1.0
 *
 * @version 1.0
 * @author Chengpu Wang
 */
public abstract class Moment {

    final public double bounding;
    final public int maxOrder;
    final public double leakage;
    final private double[] _sMoment;

    protected Moment(double bounding) {
        this.bounding = bounding;
        this.leakage = leakage(bounding);
        this._sMoment = init(bounding);
        this.maxOrder = _sMoment.length;
    }

    public double get(int n) {
        if ((n < 0) || (n >= maxOrder))
            return 0;
        return _sMoment[n];
    }

    protected abstract double[] init(double bounding);
    protected abstract double leakage(double bounding);
}
