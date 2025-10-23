import math
import os
import unittest

from fft import FFT, FFT_Order, FFT_Step, SinSource, SignalType, NoiseType, TestType
from indexSin import IndexSin

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
    fft = FFT(SinSource.Quart)

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

    def validate_Step(self, order: int, sinSource:SinSource, 
                      errPrec=5e-1, uncPrec=1e-2, precDiff=1):
        '''
        Generate dump file, and compare it with official order calculation
        Sin and Cos have significant difference for lower order
        '''
        if os.getcwd().endswith('VarianceArithemtic'):
            dumpOrderPath=f'./Python/Output/FFT_Order_{order}_{sinSource}.txt'
        elif os.getcwd().endswith('Python'):
            dumpOrderPath=f'./Output/FFT_Order_{order}_{sinSource}.txt'
        else:
            raise ValueError(f'Invalid cwd {os.getcwd()}')
        FFT_Step.dump(order, sinSource, dumpOrderPath=dumpOrderPath)
        if sinSource == SinSource.Lib:
            return
        FFT_Step.recalc(self, order, sinSource, errPrec=errPrec, uncPrec=uncPrec, precDiff=precDiff)

    def test_2(self):
        self.validate_Step(2, SinSource.Quart, errPrec=5)
        self.validate_Step(2, SinSource.Full, errPrec=5)
        self.validate_Step(2, SinSource.Prec, errPrec=5)
        self.validate_Step(2, SinSource.Lib, errPrec=5)

    def test_3(self):
        self.validate_Step(3, SinSource.Quart, errPrec=1)
        self.validate_Step(3, SinSource.Full, errPrec=1)
        self.validate_Step(3, SinSource.Prec, errPrec=1)
        self.validate_Step(3, SinSource.Lib, errPrec=1)

    def test_4(self):
        self.validate_Step(4, SinSource.Prec)
        self.validate_Step(4, SinSource.Quart,  errPrec=1)
        self.validate_Step(4, SinSource.Full)
        self.validate_Step(4, SinSource.Prec)
        self.validate_Step(4, SinSource.Lib)

    def test_5(self):
        self.validate_Step(5, SinSource.Quart)
        self.validate_Step(5, SinSource.Full)
        self.validate_Step(5, SinSource.Prec)
        self.validate_Step(5, SinSource.Lib)

    def test_6(self):
        '''
        self.validate_Step(6, SinSource.Quart)
        self.validate_Step(6, SinSource.Full)
        self.validate_Step(6, SinSource.Prec)
        '''
        self.validate_Step(6, SinSource.Lib)

    



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
            self.assertTupleEqual(tuple(ssssAggr.keys()), (SinSource.Prec, SinSource.Quart, SinSource.Full, SinSource.Lib, SinSource.Fixed))
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






if __name__ == '__main__':
    unittest.main()