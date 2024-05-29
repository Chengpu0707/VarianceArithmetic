package Type;

public class Taylor {
    
    static VarDbl[] power(double exponent) throws InitException {
        final VarDbl[] sTaylor = new VarDbl[ Momentum.MAX_FACTOR*2 ];
        sTaylor[0] = new VarDbl();
        sTaylor[1] = new VarDbl(exponent);
        --exponent;
        for (int i = 2; (i < Momentum.MAX_FACTOR*2) && (sTaylor[i - 1].value() != 0); ++i, --exponent) {
            sTaylor[i] = sTaylor[i - 1].multiply(new VarDbl(exponent / i));
        }
        return sTaylor;
    }
    
    static VarDbl[] exp() throws InitException {
        final VarDbl[] sTaylor = new VarDbl[ Momentum.MAX_FACTOR*2 ];
        sTaylor[0] = new VarDbl();
        double n = 1;
        for (int i = 1; i < Momentum.MAX_FACTOR; ++i, n *= i) {
            sTaylor[i] = new VarDbl(1.0 /n);
        }
        return sTaylor;
    }
    
    static VarDbl[] log() throws InitException {
        final VarDbl[] sTaylor = new VarDbl[ Momentum.MAX_FACTOR*2 ];
        sTaylor[0] = new VarDbl();
        for (int i = 1; i < Momentum.MAX_FACTOR*2; ++i) {
            sTaylor[i] = new VarDbl(((i%2) == 1)? +1.0/i : -1.0/i);
        }
        return sTaylor;
    }
    
    static VarDbl[] sin(double x) throws InitException {
        final VarDbl[] sTaylor = new VarDbl[ Momentum.MAX_FACTOR*2 ];
        sTaylor[0] = new VarDbl( Math.sin(x) );
        double n = 1;
        for (int i = 1; i < Momentum.MAX_FACTOR; ++i, n *= i) {
            switch (i%4) {
            case 0:
                sTaylor[i] = new VarDbl( Math.sin(x) / n ); 
                break;
            case 1:
                sTaylor[i] = new VarDbl( Math.cos(x) / n );   
                break;
            case 2:
                sTaylor[i] = new VarDbl( -Math.sin(x) / n ); 
                break;
            case 3:
                sTaylor[i] = new VarDbl( -Math.cos(x) / n ); 
                break;
            }
        }
        return sTaylor;
    }
}

