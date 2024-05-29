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
    public static final int MAX_FACTOR = 126;

    private static final VarDbl sFactor[] = new VarDbl[MAX_FACTOR];
    private static final int divid = 32;
    static {
        for (int j = 0; j < MAX_FACTOR; ++j)
            sFactor[j] = new VarDbl();
        final int limit = BINDING_FOR_TAYLOR * divid;
        final double divid2 = divid * divid;
        for (int i = -limit; i <= limit; ++i) {
            try {
                final VarDbl x2 = new VarDbl( i*i / divid2, 0);
                final VarDbl pdf = new VarDbl(1.0/Math.sqrt(2*Math.PI)).multiply(new VarDbl(Math.exp(- x2.value() * 0.5) / divid));
                for (int j = 0; j < MAX_FACTOR; ++j) {
                    sFactor[j].add(pdf, true);
                    pdf.multiply(x2, true);
                }
             } catch (InitException e) {
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
    public static VarDbl factor(int n) {
        if ((n < 0) || (n > MAX_FACTOR) || ((n % 2) == 1))
            return new VarDbl();
        return sFactor[n >> 1];
    }
}
