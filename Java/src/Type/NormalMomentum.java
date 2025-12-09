package Type;

import org.apache.commons.math3.distribution.NormalDistribution;

public class NormalMomentum extends Momentum {
    private static final NormalDistribution distr = new NormalDistribution();
 
    public NormalMomentum(double bounding) {
        super(bounding);
    }

    @Override
    protected double[] init(double bounding) {
        double term = 2 * bounding * distr.density(bounding);
        final double bounding2 = bounding * bounding;
        final int MAX_ORDER = 10000;
        final double[] sTerm = new double[MAX_ORDER];
        int n = 0;
        for (; n < MAX_ORDER; ++n) {
            sTerm[n] = term / (2*n + 1);
            if (!Double.isFinite(sTerm[n])) {
                break;
            }
            term *= bounding2;
        }
        double[] sMomentum = new double[n << 1];
        for (int i = 0; i < n; ++i) {
            sMomentum[i << 1] = sTerm[i];
        }
        for (int j = 2; j < sMomentum.length; ++j) {
            for (int i = 0; i < n; ++i) {
                sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2;
                final double prev = sMomentum[i << 1];
                sMomentum[i << 1] += sTerm[i];
                if (prev == sMomentum[i << 1]) {
                    n = i;
                    break;
                }
            }
            if (n <= 0) 
                break;
        }
        return sMomentum;
    }
    @Override
    protected double leakage(double bounding) {
        return 1 - distr.cumulativeProbability(bounding);
    }
    
}
