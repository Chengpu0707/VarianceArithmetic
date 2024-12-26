package Type;
/*
 * @since 1.0
 * 
 * @version 1.0
 * @author Chengpu Wang
 * 
 * @see IPyNb/Momentum.ipynb
 */
public class Momentum {

    public static final int BINDING_FOR_TAYLOR = 5;
    public static final int MAX_FACTOR = 442;

    private static final double sFactor[] = new double[MAX_FACTOR];
    private static final int divid = 128;
    static {
        final int limit = BINDING_FOR_TAYLOR * divid;
        final double divid2 = divid * divid;
        for (int i = -limit; i < limit; ++i) {
            try {
                final double x2 = (i + 0.5)*(i + 0.5) / divid2;
                double pdf =1.0/Math.sqrt(2*Math.PI) * Math.exp(- x2 * 0.5) / divid;
                for (int j = 0; j < MAX_FACTOR / 2; ++j) {
                    sFactor[j << 1] += pdf;
                    pdf *= x2;
                }
             } catch (Throwable e) {
                e.printStackTrace();
                assert false;
            }
        }
    }

    /*
     * Fetch variance momentum
     * 
     * @param n     order for momement factor
     * 
     * @return      0 when n is odd, variance momentum when n is even
     */
    public static double get(int n) {
        if ((n < 0) || (n >= MAX_FACTOR))
            return 0;
        return sFactor[n];
    }
}
