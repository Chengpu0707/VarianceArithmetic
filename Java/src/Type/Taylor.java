package Type;

import java.util.HashMap;
import java.util.Map;



public class Taylor {
    protected static final Map<Double, double[]> ssPowerCoeff = new HashMap<Double, double[]>();

    protected static double[] powerTaylor(double exponent) {
        final double[] sTaylor = new double[ Momentum.maxN*2 ];
        sTaylor[0] = 0;
        sTaylor[1] = exponent;
        --exponent;
        for (int i = 2; (i < Momentum.maxN*2) && (sTaylor[i - 1] != 0); ++i, --exponent) {
            sTaylor[i] = sTaylor[i - 1] * exponent / i;
        }
        return sTaylor;
    }
}

