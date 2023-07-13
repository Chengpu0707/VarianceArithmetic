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
    private static final double[] sDoubleFactorial = new double[maxFactorial*2];
    static {
        sDoubleFactorial[0] = 1;
        sDoubleFactorial[1] = 1;
        for (int i = 2; i < (maxFactorial*2); ++i) {
            sDoubleFactorial[i] = sDoubleFactorial[i-2]*i;
        }
    }
    public static double factorial(int n) {
        if ((n < 0) || (n >= maxFactorial)) return Double.NaN;
        if (n == 0) return 1;
        return sDoubleFactorial[n] * sDoubleFactorial[n - 1];
    }
    public static double doubleFactorial(int n) {
        if ((n < 0) || (n >= (maxFactorial*2))) return Double.NaN;
        return sDoubleFactorial[n];
    }

    /*
     * Calculate momentum factor
     */
    public static final int maxS = 128;
    public static final int dividS = 16;
    public static final int maxN = 200;
    private static final double ssFactor[][] = new double[maxS][maxN];
    static {
        final double factorS = 1.0 / dividS;
        for (int is = 0; is < maxS; ++is) {
            final double s = is * factorS;
            final double pdf = 2 * pdf(s);
            double power = 1;
            for (int n = 1; n <= maxN; ++n) {
                power *= s * s;
                double term = pdf * power * s / (2*n + 1);
                double sum = 0;
                double next = term;
                for (int j = 2*n + 3; sum < next; j += 2) {
                    sum = next;
                    term *= s * s / j;
                    next = sum + term;
                }
                ssFactor[is][n - 1] = sum;
            }
        }
    }
    public static double factor(int n, double s) {
        if ((n <= 0) || (maxN < n) || (s < 0)) {
            return Double.NaN;
        }
        if ((n % 2) == 1) {
            return 0; 
        }
        if ((s * dividS) >= maxS) {
            return 0;
        }
        final int lower = (int) Math.floor(s * dividS);
        final int upper = (int) Math.ceil(s * dividS);
        if (lower == upper) {
            return ssFactor[lower][n/2 - 1];
        } else if (upper < maxS) {
            return ssFactor[lower][n/2 - 1] + (ssFactor[upper][n/2 - 1] - ssFactor[lower][n/2 - 1]) * (s * dividS - lower);
        } else {
            return ssFactor[lower][n/2 - 1];
        }
    }
}
