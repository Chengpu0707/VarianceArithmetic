package Type;

import java.lang.Math;

/*
 * @since 1.0
 * 
 * @version 1.0
 * @author Chengpu Wang
 */
public class Momentum {

    public static class Normal {
        final public double bounding;
        final public int maxOrder;
        final public double leakage;
        final private double[] _sMomentum;

        // from https://introcs.cs.princeton.edu/java/22library/Gaussian.java.html
        public static double pdf(double z) {
            return Math.exp(-z*z / 2) / Math.sqrt(2 * Math.PI);
        }
        public static double cdf(double z) {
            if (z < -8.0) return 0.0;
            if (z >  8.0) return 1.0;
            double sum = 0.0, term = z;
            for (int i = 3; sum + term != sum; i += 2) {
                sum  = sum + term;
                term = term * z * z / i;
            }
            return 0.5 + sum * pdf(z);
        }

        public Normal(double bounding) {
            if ((bounding <= 0) || (10 < bounding)) {
                throw new IllegalArgumentException("Bounding must be in [1], 10]");
            }
            this.bounding = bounding;
            this.leakage = 1 - cdf(bounding);

            double term = 2 * pdf(bounding) * bounding;
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
            maxOrder = n * 2;
            _sMomentum = new double[n];
            for (int i = 0; i < n; ++i) {
                _sMomentum[i] = sTerm[i];
            }
            for (int j = 2; j < _sMomentum.length; ++j) {
                for (int i = 0; i < n; ++i) {
                    sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2;
                    final double prev = _sMomentum[i];
                    _sMomentum[i] += sTerm[i];
                    if (prev == _sMomentum[i]) {
                        n = i;
                        break;
                    }
                }
                if (n <= 0) 
                    break;
            }
        }

        public Normal() {
            this(5.0);
        }

        public double get(int n) {
            if ((n < 0) || (n >= maxOrder) || ((n & 1) != 0))
                return 0;
            return _sMomentum[n >> 1];
        }
        
    }
}
