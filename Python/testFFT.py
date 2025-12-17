import math
import os
import unittest

from fft import FFT, FFT_Signal, FFT_Order, FFT_Step, SinSource, SignalType, NoiseType, TestType
from indexSin import IndexSin
from varDbl import VarDbl

class TestFFT (unittest.TestCase):

    def testBitReverse(self):
        self.assertTupleEqual(FFT.bitReversedIndices(2), (0,2,1,3))
        self.assertTupleEqual(FFT.bitReversedIndices(3), (0,4,2,6, 1,5,3,7))
        self.assertTupleEqual(FFT.bitReversedIndices(4), (0,8,4,12, 2,10,6,14, 1,9,5,13, 3,11,7,15))

        for order in range(2, IndexSin.MAX_ORDER):
            sRes = FFT.bitReversedIndices(order)
            for i in range(len(sRes)):
                br = 0
                org = i
                for k in range(order):
                    br <<= 1
                    br |= (org & 1)
                    org >>= 1
                self.assertEqual(br, sRes[i])


class Test_FFT_Prec (unittest.TestCase):
    fft = FFT(SinSource.Prec)

    def testOrder2Sin(self):       
        sData = [0,0, 1,0, 0,0, -1,0]
        sSpec = [0,0, 0,2, 0,0, 0,-2]       
        self.assertListEqual(sSpec, Test_FFT_Prec.fft.transform(sData, True))
        self.assertListEqual(sData, Test_FFT_Prec.fft.transform(sSpec, False))

    def testOrder2Cos(self):
        sData = [1,0, 0,0, -1,0, 0,0]
        sSpec = [0,0, 2,0, 0,0, 2,0]       
        self.assertListEqual(sSpec, Test_FFT_Prec.fft.transform(sData, True))
        self.assertListEqual(sData, Test_FFT_Prec.fft.transform(sSpec, False))

    def testOrder3Sin(self):
        q2 = math.sqrt(0.5)       
        sData = [0,0, q2,0, 1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0]
        sSpec = [0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4] 
        sRes = Test_FFT_Prec.fft.transform(sData, True)
        for i, spec, res in zip(range(16), sSpec, sRes):
            try:
                self.assertAlmostEqual(spec, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  
        sRes = Test_FFT_Prec.fft.transform(sSpec, False)     
        for i, datum, res in zip(range(16), sData, sRes):
            try:
                self.assertAlmostEqual(datum, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  

    def testOrder3Cos(self):
        q2 = math.sqrt(0.5)       
        sData = [1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0, 0,0, q2,0]
        sSpec = [0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0] 
        sRes = Test_FFT_Prec.fft.transform(sData, True)
        for i, spec, res in zip(range(16), sSpec, sRes):
            try:
                self.assertAlmostEqual(spec, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  
        sRes = Test_FFT_Prec.fft.transform(sSpec, False)     
        for i, datum, res in zip(range(16), sData, sRes):
            try:
                self.assertAlmostEqual(datum, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  

    def testOrder2Linear(self):
        sData = [0,0, 1,0, 2,0, 3,0]
        sSpec = [6,0, -2,-2, -2,0, -2,2]
        sRes = Test_FFT_Prec.fft.transform(sData, True)       
        for i, spec, res in zip(range(8), sSpec, sRes):
            try:
                self.assertAlmostEqual(spec, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  
        sRes = Test_FFT_Prec.fft.transform(sSpec, False)     
        for i, datum, res in zip(range(8), sData, sRes):
            try:
                self.assertAlmostEqual(datum, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  


class Test_FFT_Quart (unittest.TestCase):
    fft = FFT(SinSource.Quart)

    def testOrder2Sin(self):       
        sData = [0,0, 1,0, 0,0, -1,0]
        sSpec = [0,0, 0,2, 0,0, 0,-2]       
        self.assertListEqual(sSpec, Test_FFT_Quart.fft.transform(sData, True))
        self.assertListEqual(sData, Test_FFT_Quart.fft.transform(sSpec, False))

    def testOrder2Cos(self):
        sData = [1,0, 0,0, -1,0, 0,0]
        sSpec = [0,0, 2,0, 0,0, 2,0]       
        self.assertListEqual(sSpec, Test_FFT_Quart.fft.transform(sData, True))
        self.assertListEqual(sData, Test_FFT_Quart.fft.transform(sSpec, False))

    def testOrder3Sin(self):
        q2 = math.sqrt(0.5)       
        sData = [0,0, q2,0, 1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0]
        sSpec = [0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4] 
        sRes = Test_FFT_Quart.fft.transform(sData, True)
        for i, spec, res in zip(range(16), sSpec, sRes):
            try:
                self.assertAlmostEqual(spec, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  
        sRes = Test_FFT_Quart.fft.transform(sSpec, False)     
        for i, datum, res in zip(range(16), sData, sRes):
            try:
                self.assertAlmostEqual(datum, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  

    def testOrder3Cos(self):
        q2 = math.sqrt(0.5)       
        sData = [1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0, 0,0, q2,0]
        sSpec = [0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0] 
        sRes = Test_FFT_Quart.fft.transform(sData, True)
        for i, spec, res in zip(range(16), sSpec, sRes):
            try:
                self.assertAlmostEqual(spec, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  
        sRes = Test_FFT_Quart.fft.transform(sSpec, False)     
        for i, datum, res in zip(range(16), sData, sRes):
            try:
                self.assertAlmostEqual(datum, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  

    def testOrder2Linear(self):
        sData = [0,0, 1,0, 2,0, 3,0]
        sSpec = [6,0, -2,-2, -2,0, -2,2]
        sRes = Test_FFT_Quart.fft.transform(sData, True)       
        for i, spec, res in zip(range(8), sSpec, sRes):
            try:
                self.assertAlmostEqual(spec, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  
        sRes = Test_FFT_Quart.fft.transform(sSpec, False)     
        for i, datum, res in zip(range(8), sData, sRes):
            try:
                self.assertAlmostEqual(datum, res, delta=math.ulp(2))
            except AssertionError as ex:
                raise ex  


@unittest.skip('No longer needed')
class Test_FFT_Full (unittest.TestCase):
    fft = FFT(SinSource.Full)

    def testOrder3Sin(self):
        q2 = math.sqrt(0.5)
        sData = [0,0, q2,0, 1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0]
        sSpec = [0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4] 
        sRes = Test_FFT_Lib.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = Test_FFT_Lib.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))

    def testOrder3Cos(self):
        q2 = math.sqrt(0.5)
        sData = [1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0, 0,0, q2,0]
        sSpec = [0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0] 
        sRes = Test_FFT_Lib.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = Test_FFT_Lib.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))


class Test_FFT_Lib (unittest.TestCase):
    fft = FFT(SinSource.Lib)

    def testOrder3Sin(self):
        q2 = math.sqrt(0.5)
        sData = [0,0, q2,0, 1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0]
        sSpec = [0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4] 
        sRes = Test_FFT_Full.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = Test_FFT_Full.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))

    def testOrder3Cos(self):
        q2 = math.sqrt(0.5)
        sData = [1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0, 0,0, q2,0]
        sSpec = [0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0] 
        sRes = Test_FFT_Full.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = Test_FFT_Full.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))


class Test_FFT_Step (unittest.TestCase):

    def test_clean(self):
        for order in range(2, 7):
            for sinSource in (SinSource.Prec, SinSource.Quart, SinSource.Lib):
                FFT_Step.dump(order, sinSource)

    def test_order3_sin1_Prec(self):
        '''
        SinSource.Prec does not mean more precision
        '''
        # i/o are accurate enough
        var = VarDbl('0.70710678118654757237', '4.83012067842292225350e-17')
        self.assertEqual(var.value(), 0.70710678118654757237)
        self.assertEqual(var.uncertainty(), 4.83012067842292225350e-17)

        fftSignal = FFT_Signal(SinSource.Prec, SignalType.Sin, 3, 1)
        fftOrder = FFT_Order(fftSignal, NoiseType.Gaussian, 0, traceSteps=True)
        self.assertEqual(fftOrder.ssSpecStep[4][10].value(), 0.0)
        self.assertEqual(fftOrder.ssSpecStep[4][10].uncertainty(), 1.18313310582018730308e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].value(), -4.44089209850062616169e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].uncertainty(), 1.18313310582018730308e-16)

        sin = fftOrder.idxSin.sin(1, 2)
        self.assertEqual(sin.value(), 7.07106781186547572737e-01)
        self.assertEqual(sin.uncertainty(), 4.83012067842292225350e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][10].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][10].uncertainty(), 6.83082217132443123058e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][11].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][11].uncertainty(), 6.83082217132443123058e-17)

        # floating point multiplication difference
        rd0 = fftOrder.ssSpecStep[3][10] * sin
        self.assertEqual(rd0.value(), 1.0000000000000002)
        self.assertEqual(rd0.uncertainty(), 8.366014421717556e-17)

    def test_order3_sin1_PrecAdj(self):
        '''
        SinSource.PrecAdj does not mean more precision
        '''
        fftSignal = FFT_Signal(SinSource.PrecAdj, SignalType.Sin, 3, 1)
        fftOrder = FFT_Order(fftSignal, NoiseType.Gaussian, 0, traceSteps=True)
        self.assertEqual(fftOrder.ssSpecStep[4][10].value(), 0.0)
        self.assertEqual(fftOrder.ssSpecStep[4][10].uncertainty(), 1.18313310582018730308e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].value(), -4.44089209850062616169e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].uncertainty(), 1.18313310582018730308e-16)

        sin = fftOrder.idxSin.sin(1, 2)
        self.assertEqual(sin.value(), 7.07106781186547572737e-01)
        self.assertEqual(sin.uncertainty(), 4.83012067842292225350e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][10].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][10].uncertainty(), 6.83082217132443123058e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][11].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][11].uncertainty(), 6.83082217132443123058e-17)

        # floating point multiplication difference
        rd0 = fftOrder.ssSpecStep[3][10] * sin
        self.assertEqual(rd0.value(), 1.0000000000000002)
        self.assertEqual(rd0.uncertainty(), 8.366014421717556e-17)

    def test_order3_sin1_Quart(self):
        '''
        similar to SinSource.Prec but with slightly higher uncertainty
        '''
        fftSignal = FFT_Signal(SinSource.Quart, SignalType.Sin, 3, 1)
        fftOrder = FFT_Order(fftSignal, NoiseType.Gaussian, 0, traceSteps=True)
        self.assertEqual(fftOrder.ssSpecStep[4][10].value(), 0.0)
        self.assertEqual(fftOrder.ssSpecStep[4][10].uncertainty(), 1.57009245868377541324e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].value(), -4.44089209850062616169e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].uncertainty(), 1.57009245868377541324e-16)

        sin = fftOrder.idxSin.sin(1, 2)
        self.assertEqual(sin.value(), 7.07106781186547572737e-01)
        self.assertEqual(sin.uncertainty(), 6.40987562127854728723e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][10].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][10].uncertainty(), 9.06493303673679024482e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][11].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][11].uncertainty(), 9.06493303673679024482e-17)

        # floating point multiplication difference
        rd0 = fftOrder.ssSpecStep[3][10] * sin
        self.assertEqual(rd0.value(), 1.0000000000000002)
        self.assertEqual(rd0.uncertainty(), 1.1102230246251568e-16)
 
        fftSignal = FFT_Signal(SinSource.Quart, SignalType.Sin, 3, 1)
        fftOrder = FFT_Order(fftSignal, NoiseType.Gaussian, 0, traceSteps=True)
        self.assertEqual(fftOrder.ssSpecStep[4][11].value(), -4.44089209850062616169e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].uncertainty(), 1.57009245868377541324e-16)

        sin = fftOrder.idxSin.sin(1, 2)
        self.assertEqual(sin.value(), 7.07106781186547572737e-01)
        self.assertEqual(sin.uncertainty(), 6.40987562127854728723e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][10].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][10].uncertainty(), 9.06493303673679024482e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][11].value(), 1.41421356237309514547e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][11].uncertainty(), 9.06493303673679024482e-17)

        # floating point multiplication difference
        rd0 = fftOrder.ssSpecStep[3][10] * sin
        self.assertEqual(rd0.value(), 1.0000000000000002)
        self.assertEqual(rd0.uncertainty(), 1.1102230246251568e-16)
 
    def test_order3_sin1_Lib(self):
        '''
        SinSource.Lib does not mean less precision
        '''
        fftSignal = FFT_Signal(SinSource.Lib, SignalType.Sin, 3, 1)
        fftOrder = FFT_Order(fftSignal, NoiseType.Gaussian, 0, traceSteps=True)
        self.assertEqual(fftOrder.ssSpecStep[4][10].value(), 2.22044604925031308085e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][10].uncertainty(), 1.57009245868377516672e-16)
        self.assertEqual(fftOrder.ssSpecStep[4][11].value(), 0.0)
        self.assertEqual(fftOrder.ssSpecStep[4][11].uncertainty(), 1.57009245868377516672e-16)

        sin = fftOrder.idxSin.sin(1, 2)
        self.assertEqual(sin.value(), 7.07106781186547572737e-01)
        self.assertEqual(sin.uncertainty(), 6.40987562127854728723e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][10].value(), 1.41421356237309492343e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][10].uncertainty(), 9.06493303673679024482e-17)
        self.assertEqual(fftOrder.ssSpecStep[3][11].value(), 1.41421356237309536752e+00)
        self.assertEqual(fftOrder.ssSpecStep[3][11].uncertainty(), 9.06493303673679024482e-17)

        # floating point multiplication difference
        rd0 = fftOrder.ssSpecStep[3][10] * sin
        self.assertEqual(rd0.value(), 1.0)
        self.assertEqual(rd0.uncertainty(), 1.1102230246251565e-16)
 



class Test_FFT_Order (unittest.TestCase):
    '''
    Check the FFT order result
    '''

    def assertErrDev(self, sOrder, forwardPrec=1e-2, reversePrec=2e-2, roundtripPrec=2e-3,
                     sNoise = [0] + [math.pow(10, n) for n in range(-17, 1)], minNoise=1e-14, minOrder=3):
        path = FFT_Order.dumpPath(sOrder=sOrder)
        if os.path.isfile(path):
            os.remove(path)
        FFT_Order.dump(sOrder, sNoise=sNoise, sNoiseType=(NoiseType.Gaussian,) )
        self.assertTrue(os.path.isfile(path))
        sssssAggr = FFT_Order.read(path)
        self.assertTrue(sssssAggr)
        self.assertTupleEqual(tuple(sssssAggr.keys()), sOrder)
        for order in sOrder:
            ssssAggr = sssssAggr[order]
            self.assertTupleEqual(tuple(ssssAggr.keys()), (SinSource.Prec, SinSource.PrecAdj, SinSource.Quart, SinSource.Lib))
            for src in ssssAggr.keys():
                sssAggr = ssssAggr[src]
                self.assertTupleEqual(tuple(sssAggr.keys()), (NoiseType.Gaussian,))
                for noiseType in (NoiseType.Gaussian,):
                    ssAggr = sssAggr[noiseType]
                    self.assertListEqual(sorted(ssAggr.keys()), sNoise)
                    for noise in sNoise:
                        sAggr = ssAggr[noise]
                        self.assertTupleEqual(tuple(sAggr.keys()), (TestType.Forward, TestType.Reverse, TestType.Roundtrip))
                        if minNoise <= noise and minOrder <= order:
                            try:
                                self.assertAlmostEqual(sAggr[TestType.Forward][0], 1, delta=forwardPrec)
                                self.assertAlmostEqual(sAggr[TestType.Reverse][0], 1, delta=reversePrec)
                                self.assertLess(0, sAggr[TestType.Roundtrip][0])
                                self.assertAlmostEqual(sAggr[TestType.Roundtrip][1], noise, delta=roundtripPrec)
                            except AssertionError as ex:
                                print(f'Error in {path} order={order} src={src} noiseType={noiseType} noise={noise}: {ex}')
                                raise ex

    def test_2to4(self):
        self.assertErrDev((2,3,4), forwardPrec=2e-1, reversePrec=5e-1, roundtripPrec=1e-2)
        
    def test_8(self):
        self.assertErrDev((8,), forwardPrec=5e-2, reversePrec=5e-2, roundtripPrec=5e-4,
                          sNoise = [0, 1e-15, 1e-12])
        
    @unittest.skip("Only for rerun")
    def test_sort(self):
        '''
        Generate filtered FFT order dump excluding specified SinSource
        Only for rerun when new SinSource is added
        '''
        FFT_Order.sort(dumpPath='./Python/Output/FFT_2_19.txt')  

    @unittest.skip("Only for rerun")
    def test_sort_filterOut(self):
        '''
        Generate filtered FFT order dump excluding specified SinSource
        Only for rerun when new SinSource is added
        '''
        FFT_Order.sort(dumpPath='./Python/Output/FFT_2_19.txt', 
                       filterFunc=lambda order, sinSource, noiseType, noise:  
                            order == 2 and sinSource == SinSource.Lib)






if __name__ == '__main__':
    unittest.main()