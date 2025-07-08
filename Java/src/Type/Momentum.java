package Type;


/*
 * @since 1.0
 * 
 * @version 1.0
 * @author Chengpu Wang
 */
public abstract class Momentum {
    
    final public double bounding;
    final public int maxOrder;
    final public double leakage;
    final private double[] _sMomentum;

    protected Momentum(double bounding) {
        this.bounding = bounding;
        this.leakage = leakage(bounding);
        this._sMomentum = init(bounding);
        this.maxOrder = _sMomentum.length;
    }

    public double get(int n) {
        if ((n < 0) || (n >= maxOrder))
            return 0;
        return _sMomentum[n];
    }

    protected abstract double[] init(double bounding);
    protected abstract double leakage(double bounding);
}


