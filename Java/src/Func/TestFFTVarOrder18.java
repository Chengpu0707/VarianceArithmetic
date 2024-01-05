package Func;

import org.junit.Test;

public class TestFFTVarOrder18 {
    static final RealType realType = RealType.Var;
    static final int order = 18;
    
    @Test
    public void testSinFreq1() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 1);
    }


    @Test
    public void testSinFreq2() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 2);
    }

    
}
