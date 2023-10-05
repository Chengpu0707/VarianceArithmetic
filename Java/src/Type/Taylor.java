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
    
    static double[] log() {
        final double[] sTaylor = new double[ Momentum.maxN*2 ];
        sTaylor[0] = 0;
        for (int i = 1; i < Momentum.maxN*2; ++i) {
            sTaylor[i] = ((i%2) == 1)? +1.0/i : -1.0/i;
        }
        return sTaylor;
    }
    
    static double[] sin(double x) {
        final double[] sTaylor = new double[ Momentum.maxN*2 ];
        sTaylor[0] = Math.sin(x);
        double n = 1;
        for (int i = 1; i < Momentum.maxFactorial; ++i, n *= i) {
            switch (i%4) {
            case 0:
                sTaylor[i] = Math.sin(x) / n; 
                break;
            case 1:
                sTaylor[i] = Math.cos(x) / n;   
                break;
            case 2:
                sTaylor[i] = -Math.sin(x) / n; 
                break;
            case 3:
                sTaylor[i] = -Math.cos(x) / n; 
                break;
            }
        }
        return sTaylor;
    }
}

