package Func;

import org.junit.Test;

public class TestFFTVarOrder4 {
    static final RealType realType = RealType.Var;
    static final int order = 4;
    
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
    public void testSinFreq5() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Sin, order, 5);
    }

    @Test
    public void testSinFreq6() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
              SignalType.Sin, order, 6);
    }

    @Test
    public void testSinFreq7() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
              SignalType.Sin, order, 7);
    }

    @Test
    public void testCosFreq1() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 1);
    }


    @Test
    public void testCosFreq2() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 2);
    }

    @Test
    public void testCosFreq3() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 3);
    }
    
    @Test
    public void testCosFreq4() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 4);
    }


    @Test
    public void testCosFreq5() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 5);
    }

    @Test
    public void testCosFreq6() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 6);
    }
    
    @Test
    public void testCosFreq7() {
        TestFFT.test(realType, NoiseType.Gaussian, 0, 
                SignalType.Cos, order, 7);
    }
    
}
