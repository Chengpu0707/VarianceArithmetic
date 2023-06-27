package Type;

/*
 * copy from https://introcs.cs.princeton.edu/java/22library/Gaussian.java.html#:~:text=PI%29%3B%7D%2F%2F%20return%20pdf%28x%2C%20mu%2C%20sigma%29%20%3D%20Gaussian%20pdf,bisection%20search%40DeprecatedpublicstaticdoublePhiInverse%28doubley%29%7BreturninverseCDF%28y%29%3B%7D%2F%2F%20test%20clientpublicstaticvoidmain%28String%5B%5Dargs%29%7Bdoublez%20%3DDouble.parseDouble%28args%29%3Bdoublemu%20%3DDouble.parseDouble%28args%29%3Bdoublesigma%20%3DDouble.parseDouble%28args%29%3BStdOut.println%28cdf%28z%2Cmu%2Csigma%29%29%3Bdoubley%20%3Dcdf%28z%29%3BStdOut.println%28inverseCDF%28y%29%29%3B%7D%7D
 * 
 * 
 */
public class Momentum {

    // return pdf(x) = standard Gaussian pdf
    public static double pdf(double z) {
        return Math.exp(-z*z / 2) / Math.sqrt(2 * Math.PI);
    }
    
    // return cdf(z) = standard Gaussian cdf using Taylor approximation
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

    // Compute z such that cdf(z) = y via bisection search
    public static double inverseCDF(double y) {
        return inverseCDF(y, 0.00000001, -8, 8);
    }

    // bisection search
    private static double inverseCDF(double y, double delta, double lo, double hi) {
        double mid = lo + (hi - lo) / 2;
        if (hi - lo < delta) return mid;
        if (cdf(mid) > y) return inverseCDF(y, delta, lo, mid);
        else              return inverseCDF(y, delta, mid, hi);
    }

    public static final int maxFactorial = 171;
    private static final double[] sDoubleFactorial = new double[maxFactorial];
    static {
        sDoubleFactorial[0] = 1;
        sDoubleFactorial[1] = 1;
        for (int i = 2; i < maxFactorial; ++i) {
            sDoubleFactorial[i] = sDoubleFactorial[i-2]*i;
        }
    }
    public static double factorial(int n) {
        if ((n < 0) || (n >= maxFactorial)) return Double.NaN;
        if (n == 0) return 1;
        return sDoubleFactorial[n] * sDoubleFactorial[n - 1];
    }
    public static double doubleFactorial(int n) {
        if ((n < 0) || (n >= maxFactorial)) return Double.NaN;
        return sDoubleFactorial[n];
    }

    public static final int maxX = 512;
    public static final int dividX = 16;
    public static final int maxN = 100;    // maxFactorial = 171
    private static final double ssResidual[][] = new double[maxX][maxN];
    static {
        final double factorX = 1.0 / dividX;
        for (int ix = 0; ix < maxX; ++ix) {
            double x = ix * factorX;
            for (int n = 1; n <= maxN; ++n) {
                final boolean odd = ((n % 2) == 1);
                double term = odd?  1 : x;
                double sum = term;
                for (int j = odd? 2 : 3; j < n; j += 2) {
                    term *= x * x / j;
                    sum += term;
                }
                ssResidual[ix][n - 1] = sum * pdf(x);
            }
        }
    }
    public static double residual(int n, double z) {
        if ((n <= 0) || (maxN < n) || (z < 0)) {
            return Double.NaN;
        }
        if ((z * dividX) >= maxX) {
            return 0;
        }
        final int lower = (int) Math.floor(z * dividX);
        final int upper = (int) Math.ceil(z * dividX);
        if (upper < maxX) {
            return ssResidual[lower][n - 1] + (ssResidual[upper][n - 1] - ssResidual[lower][n - 1]) * (z * dividX - lower);
        } else {
            return ssResidual[lower][n - 1];
        }
    }

    public static double factor(int n, double z) {
        if (Double.isFinite(z)) {
            return ((n % 2) == 1)? residual(n, z) : cdf(z) - residual(n, z);
        } else {
            return ((n % 2) == 1)? 0 : 1;
        }
    }
}
