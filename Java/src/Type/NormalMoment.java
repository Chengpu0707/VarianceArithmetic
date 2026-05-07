package Type;

/**
 * Normal-distribution implementation of {@link Moment}: computes the truncated
 * Gaussian moment table up to the given bounding factor for use in
 * statistical Taylor expansion.
 */
import org.apache.commons.math3.distribution.NormalDistribution;

public class NormalMoment extends Moment {
    private static final NormalDistribution distr = new NormalDistribution();

    public NormalMoment(double bounding) {
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
        double[] sMoment = new double[n << 1];
        for (int i = 0; i < n; ++i) {
            sMoment[i << 1] = sTerm[i];
        }
        for (int j = 2; j < sMoment.length; ++j) {
            for (int i = 0; i < n; ++i) {
                sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2;
                final double prev = sMoment[i << 1];
                sMoment[i << 1] += sTerm[i];
                if (prev == sMoment[i << 1]) {
                    n = i;
                    break;
                }
            }
            if (n <= 0)
                break;
        }
        // Normalize per Formula (2.2): divide by ∫ρ dz = 1 - leakage so ζ(0,κ)=1.
        final double normFactor = 1 - leakage(bounding);
        for (int i = 0; i < sMoment.length; ++i) {
            sMoment[i] /= normFactor;
        }
        return sMoment;
    }
    @Override
    protected double leakage(double bounding) {
        // Two-sided leakage = 1 - ∫_{-κ}^{κ} ρ(z) dz = 2 · P(Z > κ).
        return 2 * (1 - distr.cumulativeProbability(bounding));
    }

}
