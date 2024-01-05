package Func;

import org.junit.Test;

public class TestFFTVarSlope {
    static final RealType realType = RealType.Var;
    
    @Test
    public void testSlope2() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Slope, 2, 0);
    }

    @Test
    public void testSlope3() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Slope, 3, 0);
    }

    @Test
    public void testSlope4() {
        TestFFT.test(realType, NoiseType.Gaussian, 1E-2, 
                SignalType.Slope, 4, 0);
    }

    @Test
    public void testSlope5() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Slope, 5, 0);
    }
}
