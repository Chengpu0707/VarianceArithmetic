package Type;




public class Taylor {
    
    static double[] power(double exponent) {
        final double[] sTaylor = new double[ Momentum.maxN*2 ];
        sTaylor[0] = 0;
        sTaylor[1] = exponent;
        --exponent;
        for (int i = 2; (i < Momentum.maxN*2) && (sTaylor[i - 1] != 0); ++i, --exponent) {
            sTaylor[i] = sTaylor[i - 1] * exponent / i;
        }
        return sTaylor;
    }
    
    static double[] exp() {
        final double[] sTaylor = new double[ Momentum.maxN*2 ];
        sTaylor[0] = 0;
        double n = 1;
        for (int i = 1; i < Momentum.maxFactorial; ++i, n *= i) {
            sTaylor[i] = 1.0/n;
        }
        return sTaylor;
    }
}

