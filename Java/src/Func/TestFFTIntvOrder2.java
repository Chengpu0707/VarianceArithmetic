package Func;

import org.junit.Test;

public class TestFFTIntvOrder2 {
    static final RealType realType = RealType.Intv;
    static final int order = 2;

    @Test
    public void testSinFreq1() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 1);
    }

    @Test
    public void testCosFreq1() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 1);
    }
}
