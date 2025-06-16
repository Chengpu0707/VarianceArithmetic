import math
import unittest

from fft import FFTBase, FFTIndexSin, FFTLibSin, FFTTest

class TestFFTBase (unittest.TestCase):

    def testBitReverse(self):
        self.assertTupleEqual(FFTBase.bitReversedIndices(2), (0,2,1,3))
        self.assertTupleEqual(FFTBase.bitReversedIndices(3), (0,4,2,6, 1,5,3,7))
        self.assertTupleEqual(FFTBase.bitReversedIndices(4), (0,8,4,12, 2,10,6,14, 1,9,5,13, 3,11,7,15))

        for order in range(2, FFTBase.MAX_ORDER):
            sRes = FFTBase.bitReversedIndices(order)
            for i in range(len(sRes)):
                br = 0
                org = i
                for k in range(order):
                    br <<= 1
                    br |= (org & 1)
                    org >>= 1
                self.assertEqual(br, sRes[i])


class TestFFTIndexSin (unittest.TestCase):
    fft = FFTIndexSin()

    def testOrder2Sin(self):
        self.assertEqual(TestFFTIndexSin.fft.sin(-1,2), -1)
        self.assertEqual(TestFFTIndexSin.fft.sin(0,2), 0)
        self.assertEqual(TestFFTIndexSin.fft.sin(1,2), 1)
        self.assertEqual(TestFFTIndexSin.fft.sin(2,2), 0)
        self.assertEqual(TestFFTIndexSin.fft.sin(3,2), -1)
        self.assertEqual(TestFFTIndexSin.fft.sin(4,2), 0)
        
        sData = [0,0, 1,0, 0,0, -1,0]
        sSpec = [0,0, 0,2, 0,0, 0,-2]       
        self.assertListEqual(sSpec, TestFFTIndexSin.fft.transform(sData, True))
        self.assertListEqual(sData, TestFFTIndexSin.fft.transform(sSpec, False))

    def testOrder2Cos(self):
        self.assertEqual(TestFFTIndexSin.fft.cos(-1,2), 0)
        self.assertEqual(TestFFTIndexSin.fft.cos(0,2), 1)
        self.assertEqual(TestFFTIndexSin.fft.cos(1,2), 0)
        self.assertEqual(TestFFTIndexSin.fft.cos(2,2), -1)
        self.assertEqual(TestFFTIndexSin.fft.cos(3,2), 0)
        self.assertEqual(TestFFTIndexSin.fft.cos(4,2), 1)

        sData = [1,0, 0,0, -1,0, 0,0]
        sSpec = [0,0, 2,0, 0,0, 2,0]       
        self.assertListEqual(sSpec, TestFFTIndexSin.fft.transform(sData, True))
        self.assertListEqual(sData, TestFFTIndexSin.fft.transform(sSpec, False))

    def testOrder3Sin(self):
        q2 = math.sqrt(0.5)
        self.assertEqual(TestFFTIndexSin.fft.sin(-1,3), -q2)
        self.assertEqual(TestFFTIndexSin.fft.sin(0,3), 0)
        self.assertEqual(TestFFTIndexSin.fft.sin(1,3), q2)
        self.assertEqual(TestFFTIndexSin.fft.sin(2,3), 1)
        self.assertEqual(TestFFTIndexSin.fft.sin(3,3),  q2)
        self.assertEqual(TestFFTIndexSin.fft.sin(4,3), 0)
        self.assertEqual(TestFFTIndexSin.fft.sin(5,3), -q2)
        self.assertEqual(TestFFTIndexSin.fft.sin(6,3), -1)
        self.assertEqual(TestFFTIndexSin.fft.sin(7,3), -q2)
        self.assertEqual(TestFFTIndexSin.fft.sin(8,3), 0)
        
        sData = [0,0, q2,0, 1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0]
        sSpec = [0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4] 
        sRes = TestFFTIndexSin.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = TestFFTIndexSin.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))

    def testOrder3Cos(self):
        q2 = math.sqrt(0.5)
        self.assertEqual(TestFFTIndexSin.fft.cos(-1,3), q2)
        self.assertEqual(TestFFTIndexSin.fft.cos(0,3), 1)
        self.assertEqual(TestFFTIndexSin.fft.cos(1,3), q2)
        self.assertEqual(TestFFTIndexSin.fft.cos(2,3), 0)
        self.assertEqual(TestFFTIndexSin.fft.cos(3,3), -q2)
        self.assertEqual(TestFFTIndexSin.fft.cos(4,3), -1)
        self.assertEqual(TestFFTIndexSin.fft.cos(5,3), -q2)
        self.assertEqual(TestFFTIndexSin.fft.cos(6,3), 0)
        self.assertEqual(TestFFTIndexSin.fft.cos(7,3), q2)
        self.assertEqual(TestFFTIndexSin.fft.cos(8,3), 1)
        
        sData = [1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0, 0,0, q2,0]
        sSpec = [0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0] 
        sRes = TestFFTIndexSin.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = TestFFTIndexSin.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))


class TestFFTLibSin (unittest.TestCase):
    fft = FFTLibSin()

    def testOrder3Sin(self):
        q2 = math.sqrt(0.5)
        error = 1.2246467991473532e-16
        self.assertEqual(TestFFTLibSin.fft.sin(-1,3), -q2)
        self.assertEqual(TestFFTLibSin.fft.sin(0,3), 0)
        self.assertEqual(TestFFTLibSin.fft.sin(1,3), q2)
        self.assertEqual(TestFFTLibSin.fft.sin(2,3), 1)
        self.assertEqual(TestFFTLibSin.fft.sin(3,3),  q2)
        self.assertEqual(TestFFTLibSin.fft.sin(4,3), error)
        self.assertEqual(TestFFTLibSin.fft.sin(5,3), -q2 + error)
        self.assertEqual(TestFFTLibSin.fft.sin(6,3), -1)
        self.assertEqual(TestFFTLibSin.fft.sin(7,3), -q2 - error)
        self.assertEqual(TestFFTLibSin.fft.sin(8,3), -error * 2)
        
        sData = [0,0, q2,0, 1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0]
        sSpec = [0,0, 0,4, 0,0, 0,0, 0,0, 0,0, 0,0, 0,-4] 
        sRes = TestFFTLibSin.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = TestFFTLibSin.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))

    def testOrder3Cos(self):
        q2 = math.sqrt(0.5)
        error = 1.2246467991473532e-16
        self.assertEqual(TestFFTLibSin.fft.cos(-1,3), q2)
        self.assertEqual(TestFFTLibSin.fft.cos(0,3), 1)
        self.assertEqual(TestFFTLibSin.fft.cos(1,3), q2)
        self.assertEqual(TestFFTLibSin.fft.cos(2,3), error/2)
        self.assertEqual(TestFFTLibSin.fft.cos(3,3), -q2 + error/2)
        self.assertEqual(TestFFTLibSin.fft.cos(4,3), -1)
        self.assertEqual(TestFFTLibSin.fft.cos(5,3), -q2 - math.ulp(1)/2)
        self.assertEqual(TestFFTLibSin.fft.cos(6,3), -error *3/2)
        self.assertEqual(TestFFTLibSin.fft.cos(7,3), q2 - math.ulp(1))
        self.assertEqual(TestFFTLibSin.fft.cos(8,3), 1)
        
        sData = [1,0, q2,0, 0,0, -q2,0, -1,0, -q2,0, 0,0, q2,0]
        sSpec = [0,0, 4,0, 0,0, 0,0, 0,0, 0,0, 0,0, 4,0] 
        sRes = TestFFTLibSin.fft.transform(sData, True)
        for spec, res in zip(sSpec, sRes):
            self.assertAlmostEqual(spec, res, delta=math.ulp(2))     
        sRes = TestFFTLibSin.fft.transform(sSpec, False)     
        for datum, res in zip(sData, sRes):
            self.assertAlmostEqual(datum, res, delta=math.ulp(2))


class TestDumpFFT (unittest.TestCase):

    def testOrder_4(self):
        FFTTest.dumpOrders(sOrder=[4])

    @unittest.skip('Too slow')
    def testOrderAll(self):
        FFTTest.dumpOrders()

    def testSpectra(self):
        '''
        Demonstrate the reverse has not enough calculation for large enough uncertainty
        '''
        with open(f'./Python/Output/FFT_4_8_Spec.txt', 'w') as fw:
            FFTTest.dumpSpectrumHeader(fw)
            FFTTest.dumpSpectra(fw, range(4,8))




if __name__ == '__main__':
    unittest.main()