package Type;

public class Momentum {

    // return pdf(x) = standard Gaussian pdf
    public static double pdf(double z) {
        return Math.exp(-z*z / 2) / Math.sqrt(2 * Math.PI);
    }
    
    public static final int maxFactorial = 200;
    private static final double[] sDoubleFactorial = new double[maxFactorial*2];
    static {
        sDoubleFactorial[0] = 1;
        sDoubleFactorial[1] = 1;
        for (int i = 2; i < (maxFactorial*2); ++i) {
            sDoubleFactorial[i] = sDoubleFactorial[i-2]*i;
        }
    }
    public static double doubleFactorial(int n) {
        if ((n < 0) || (n >= (maxFactorial*2))) return Double.NaN;
        return sDoubleFactorial[n];
    }

    /*
     * Calculate momentum factor
     */
    public static final int BINDING_FOR_TAYLOR = 5;

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
    public static double factor(int n, double s, boolean normalized) {
        if (n == 0)
            return normalized? 1.0 : 0.99999947194737793;
        if ((n < 0) || (maxN*2 < n) || (s < 0)) {
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
        final double ret;
        if (lower == upper) {
            ret = ssFactor[lower][n/2 - 1];
        } else if (upper < maxS) {
            ret = ssFactor[lower][n/2 - 1] + (ssFactor[upper][n/2 - 1] - ssFactor[lower][n/2 - 1]) * (s * dividS - lower);
        } else {
            ret = ssFactor[lower][n/2 - 1];
        }
        if (normalized) {
            if (n == 2)
                return 1;
            else
                return ret / factor(2, s, false); 
        } else
            return ret; 
    }
    public static double factor(int n, double s) {
        return factor(n, s, true);
    }
    public static double factor(int n) {
        return factor(n, BINDING_FOR_TAYLOR, false);
    }
}
