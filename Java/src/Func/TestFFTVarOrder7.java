package Func;

import org.junit.Test;

public class TestFFTVarOrder7 {
    static final RealType realType = RealType.Var;
    static final int order = 7;
    
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
    
    @Test
    public void testSinFreq3() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 3);
    }
    
    @Test
    public void testSinFreq4() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 4);
    }
    
    @Test
    public void testSinFreq15() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 15);
    }
    
    @Test
    public void testSinFreq15_noisy15() {
        TestFFT.test(realType, NoiseType.Gaussian, 1e-15, 
                SignalType.Sin, order, 15);
    }
    
    @Test
    public void testSinFreq15_noisy16() {
        TestFFT.test(realType, NoiseType.Gaussian, 4e-16, 
                SignalType.Sin, order, 15);
    }
    
    @Test
    public void testSinFreq31() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 31);
    }
    
}
