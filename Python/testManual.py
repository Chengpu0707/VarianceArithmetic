'''
Unit test that needs more than one minutes to finish.

To run a particular test:
    VarianceArithmetic\Python> python -m unittest testManual.<class>.<method>
'''
from collections.abc import Callable
import datetime
import functools
import logging
import math
import numpy
import os
import random
import os
import scipy
import time
import unittest

from fft import NoiseType, FFT_Order
from indexSin import OUTDIR
from matrix import createHilbertMatrix, addNoise
from matrix import adjugate
import momentum
from regressiveSin import RegressiveSin
from taylor import Taylor, NotMonotonicException, NotStableException
import taylor
from varDbl import VarDbl
import varDbl

from testMatrix import Adjugate, TestAdjugate


SKIP_TEST = True

logger = logging.getLogger(__name__)


MAX_COUNT = 100000

def calcStat(sSample:numpy.array, samples:int) -> tuple[tuple[float, float]]:
    '''
    Reuse the random sampling array {sSample} to calculate the sample mean and standard deviation
    '''
    return tuple([(sSample[i: i+samples].mean(), sSample[i: i+samples].std())
                  for i in range(0, min(len(sSample), MAX_COUNT*samples), samples)])


class TestNormal (unittest.TestCase):
    INV_SQRT2 = 1.0 / math.sqrt(2)
    SAMPLES = (2,3,5,7,10,20,30,50,70,100,200,500,1000,2000,5000,10000,20000,50000)
    KAPPAS = (1,1.25,1.5,1.75, 2,2.25,2.5,2.75, 3,3.5, 4,4.5, 5,5.5, 6)
    MIN_COUNT = 10000

    header = 'Samples\tKappa\tCount\tMean\tDeviation\tRange\n'
    leakPath = f"{OUTDIR}/Python/Output/NormalLeakage.txt"
    resultPath = f"{OUTDIR}/Python/Output/NormalSamples.txt"

    sSample = None

    @staticmethod
    def initSamples():
        if not TestNormal.sSample:
            print(f'Start generating {TestNormal.MIN_COUNT} * {TestNormal.SAMPLES[-1]} samples for normal distribution at {datetime.datetime.now()}')
            TestNormal.sSample = scipy.stats.norm.rvs(size = TestNormal.MIN_COUNT * TestNormal.SAMPLES[-1])


    @staticmethod
    def inputLeak(mu:float, sigma:float, kapp:float):
        '''
        The sample bounding leakage of normal distribution for the bounding factor {kappa}.
        {mu} is the sample mean, {sigma} is the sample standard deviation.
        '''
        return 1 - scipy.special.erf(abs(kapp*sigma - mu) * TestNormal.INV_SQRT2)*0.5 \
                 - scipy.special.erf(abs(kapp*sigma + mu) * TestNormal.INV_SQRT2)*0.5

    def test_inputLeak(self):
        for k in range(1, 7):
            self.assertAlmostEqual(k, scipy.special.erfinv(1 - TestNormal.inputLeak(0, 1, k)) * math.sqrt(2))
            
    #@unittest.skipIf(SKIP_TEST, 'Ran 1 tests in 83s')
    def test_input(self):
        TestNormal.initSamples()
        sPrev = {}
        with open(TestNormal.leakPath, 'w') as f:
            f.write(TestNormal.header)
            prev = 1
            for k in TestNormal.KAPPAS:
                sPrev[k] = TestNormal.inputLeak(0, 1, k)
                self.assertLess(sPrev[k], prev)
                prev = sPrev[k]
                f.write(f'0\t{k}\t0\t{prev}\t0\t{k}\n')
            f.flush()
            for samples in reversed(TestNormal.SAMPLES):
                print(f'Start calculate samples={samples} for Normal distribution at {datetime.datetime.now()}')
                sStat = calcStat(TestNormal.sSample, samples)
                prev = 1
                for k in TestNormal.KAPPAS:
                    sLeak = [TestNormal.inputLeak(mu, sigma, k) for mu, sigma in sStat]
                    mean = numpy.mean(sLeak)
                    self.assertGreater(mean, sPrev[k])
                    self.assertLess(mean, prev)
                    prev = sPrev[k] = mean
                    kapp = scipy.special.erfinv(1 - mean) * math.sqrt(2)
                    self.assertLess(kapp, k)
                    f.write(f'{samples}\t{k}\t{len(sLeak)}\t{mean}\t{numpy.std(sLeak)}\t{kapp}\n')
                    f.flush()
        print(f'Finsh writing to {TestNormal.leakPath}')

    @staticmethod
    def outputLeak(mu:float, sigma:float, kapp:float):
        '''
        The result leakage of normal distribution for the bounding factor {kappa}.
        {mu} is the sample mean, {sigma} is the sample standard deviation.
        '''
        k1 = abs(mu + kapp*sigma)
        k2 = abs(mu - kapp*sigma)
        return 1 - scipy.special.erf(k1 * TestNormal.INV_SQRT2)*0.5 + k1*scipy.stats.norm.pdf(k1) \
                 - scipy.special.erf(k2 * TestNormal.INV_SQRT2)*0.5 + k2*scipy.stats.norm.pdf(k2)
    
    def test_outputLeak(self):
        self.assertAlmostEqual(TestNormal.outputLeak(0, 1, 5), 1.544050e-05)
        for k in TestNormal.KAPPAS:
            try:
                self.assertLess(TestNormal.outputLeak(0, 1, k), TestNormal.outputLeak(0.1, 1, k))
                if k < 1.5:
                    self.fail(f'kappa={k} should not pass the test')
            except AssertionError:
                if k < 1.5:
                    continue
                raise
     
    @unittest.skipIf(SKIP_TEST, 'Ran 1 test in 3171s')
    def test_output(self):
        '''
        extract the sample count from 1/kappa.
        '''
        TestNormal.initSamples()
        RESOLUTION = 1000
        MIN_KAPPA = 1.5
        sPrev = {}
        with open(TestNormal.resultPath, 'w') as f:
            f.write(TestNormal.header)
            prev = 1
            for k in TestNormal.KAPPAS:
                sPrev[k] = TestNormal.outputLeak(0, 1, k)
                self.assertGreater(sPrev[k], 0)
                self.assertLess(sPrev[k], 1)
                self.assertLess(sPrev[k], prev)
                prev = sPrev[k]
                f.write(f'0\t{k}\t0\t{prev}\t0\t{k}\n')
            f.flush()
            for samples in reversed(TestNormal.SAMPLES):
                print(f'Start calculate samples={samples} for Normal distribution at {datetime.datetime.now()}')
                sStat = calcStat(TestNormal.sSample, samples)
                prev = 1
                for k in TestNormal.KAPPAS:
                    if k <= MIN_KAPPA:
                        continue
                    sLeak = [TestNormal.outputLeak(mu, sigma, k) for mu, sigma in sStat]
                    mean = numpy.mean(sLeak)
                    self.assertGreater(mean, 0)
                    self.assertLess(mean, 1)
                    self.assertGreater(mean, sPrev[k])
                    self.assertLess(mean, prev)
                    prev = sPrev[k] = mean
                    left = leftest = int(MIN_KAPPA * RESOLUTION)
                    right = rightest = int(k * RESOLUTION)
                    while (left + 1) < right:
                        mid = (left + right) // 2
                        testLeak = TestNormal.outputLeak(0, 1, mid/RESOLUTION)
                        if testLeak > mean:
                            left = mid
                        else:
                            right = mid
                    if left == leftest:
                        continue
                    if right == rightest:
                        continue
                    kappa = left / RESOLUTION
                    right /= RESOLUTION
                    self.assertGreater(k, kappa)
                    f.write(f'{samples}\t{k}\t{len(sLeak)}\t{mean}\t{numpy.std(sLeak)}\t{kappa}\n')
                    f.flush()
        print(f'Finsh writing to {TestNormal.resultPath}')


class TestUniform (unittest.TestCase):
    SQRT3 = math.sqrt(3)
    SAMPLES = (2,3,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)
    MIN_COUNT = 10000

    header = 'Samples\tCount\tMean\tDeviation\n'
    inputPath = f'{OUTDIR}/Python/Output/UniformLeakage.txt'
    outputPath = f'{OUTDIR}/Python/Output/UniformSamples.txt'

    sSample = None

    @staticmethod
    def initSamples():
        if not TestUniform.sSample:
            print(f'Start generating {TestUniform.MIN_COUNT} * {TestUniform.SAMPLES[-1]} samples for uniform distribution at {datetime.datetime.now()}')
            TestUniform.sSample = scipy.stats.uniform.rvs(loc = -TestUniform.SQRT3, scale = 2*TestUniform.SQRT3, 
                                                          size = TestUniform.MIN_COUNT * TestUniform.SAMPLES[-1])   

    @staticmethod
    def inputLeak(mu:float, sigma:float):
        '''
        The sample bounding leakage of uniform distribution between [0,1].
        {mu} is the sample mean, {sigma} is the sample standard deviation.
        '''
        sigma *= TestUniform.SQRT3
        err = 0.0
        for k in (abs(mu + sigma), abs(mu - sigma)):
            if k > TestUniform.SQRT3:
                err += (k - TestUniform.SQRT3) 
        return err / (2*TestUniform.SQRT3)

    def test_uniformLeak(self):
        self.assertAlmostEqual(TestUniform.inputLeak(0, 1), 0)
        self.assertAlmostEqual(TestUniform.inputLeak(0.2 * TestUniform.SQRT3, 1), 0.1)
        self.assertAlmostEqual(TestUniform.inputLeak(-0.2 * TestUniform.SQRT3, 1), 0.1)
        self.assertAlmostEqual(TestUniform.inputLeak(0, 1.1), 0.1)
        self.assertAlmostEqual(TestUniform.inputLeak(0, 0.9), 0)
        self.assertAlmostEqual(TestUniform.inputLeak(0.1 * TestUniform.SQRT3, 1.1), 0.1)

    @staticmethod
    def outputLeak(mu:float, sigma:float):
        '''
        The sample bounding leakage of uniform distribution between [0,1].
        {mu} is the sample mean, {sigma} is the sample standard deviation.
        '''
        return abs(sigma**3 - 1)

    @unittest.skipIf(SKIP_TEST, 'Ran 1 tests in 75s')
    def test(self):
        TestUniform.initSamples()
        with open(TestUniform.inputPath, 'w') as f0, open(TestUniform.outputPath, 'w') as f1:
            f0.write(TestUniform.header)
            f1.write(TestUniform.header)
            sPrev = [0, 0]
            for samples in reversed(TestUniform.SAMPLES):
                print(f'Start calculate samples={samples} for uniform distribution at {datetime.datetime.now()}')
                sStat = calcStat(TestUniform.sSample, samples)
                ssLeak = [[TestUniform.inputLeak(mu, sigma) for mu, sigma in sStat],
                          [TestUniform.outputLeak(mu, sigma) for mu, sigma in sStat]]
                for i, sLeak in enumerate(ssLeak):
                    mean = numpy.mean(sLeak)
                    if i == 0:
                        self.assertGreater(mean, 0)
                        self.assertLess(mean, 1)
                        self.assertGreater(mean, sPrev[i])
                        sPrev[i] = mean
                        f0.write(f'{samples}\t{len(sLeak)}\t{mean}\t{numpy.std(sLeak)}\n')
                        f0.flush()
                    else:
                        f1.write(f'{samples}\t{len(sLeak)}\t{abs(mean)}\t{numpy.std(sLeak)}\n')
                        f1.flush()
 

class TestBoundingRange (unittest.TestCase):

    def assert_func(self, fileName:str, sKappa=(2,2.5, 3,3.5, 4,4.5, 5,5.5, 6)):
        sMomentum = {}
        ssBounding = {}
        filePath = f'{OUTDIR}/Python/Output/{fileName}Samples.txt'
        match fileName:
            case 'Normal':
                with open(filePath) as f:
                    hdr = next(f)
                    self.assertEqual(hdr, 'Samples\tKappa\tCount\tMean\tDeviation\tRange\n')
                    for line in f:
                        n, k, cnt, leak, std, bounding = map(float, line.split('\t'))
                        if k not in sKappa:
                            continue
                        sBounding = ssBounding.setdefault(k, {})
                        sBounding[int(n)] = bounding
                        if bounding not in sMomentum:
                            sMomentum[bounding] = momentum.Normal(bounding=bounding)
            case _:
                raise ValueError(f'Invalid fileName {fileName}')

        sFunc = {
            'x': lambda mmt, var: taylor.Taylor.polynominal1d(var, (0,1), momentum=mmt),
            'sin(x)': lambda mmt, var: taylor.Taylor.sin(var, momentum=mmt),
            'exp(x)': lambda mmt, var: taylor.Taylor.exp(var, momentum=mmt), 
            'log(x)': lambda mmt, var: taylor.Taylor.log(var, momentum=mmt), 
        }
        var = varDbl.VarDbl(1, 0.1)
        zero = varDbl.VarDbl(1, 0)

        def powFunc(exp, mmt, var):
            return taylor.Taylor.pow(var, exp, momentum=mmt)

        for exp in (2, 0.5, -1, -2):
            sFunc[f'x^{exp}'] = functools.partial(powFunc, exp)

        filePath = f'{OUTDIR}/Python/Output/{fileName}Bounding.txt'
        with open(filePath, 'w') as f:
            f.write('Kappa\tSamples\tBounding\tLeakage\tFunction\tStable Variance\tOutput Variance\tVariance Ratio\tVariance Leak\n')
            for calc, func in sFunc.items():
                stable = func(sMomentum[6], var)
                for kappa_s, sBounding in ssBounding.items():
                    for samples, bounding in sBounding.items():
                        mmt = sMomentum[bounding]
                        res = func(mmt, var)
                        f.write(f'{kappa_s}\t{samples}\t{bounding}\t{mmt.leakage}\t{calc}')
                        f.write(f'\t{stable.variance()}\t{res.variance()}')
                        f.write(f'\t{res.variance()/stable.variance()}\t{1 - res.variance()/stable.variance()}\n')

    @unittest.skipIf(SKIP_TEST, 'Ran 1 tests in 178s')
    def test_Normal(self):
        self.assert_func('Normal')


class Test_FFT_Order (unittest.TestCase):

    @unittest.skipIf(SKIP_TEST, 'Too slow')
    def test_Gaussian(self):
        FFT_Order.dump(sNoiseType=(NoiseType.Gaussian,))

    @unittest.skipIf(SKIP_TEST, 'Too slow')
    def test_White(self):
        FFT_Order.dump(sNoiseType=(NoiseType.White,))



class TestAdjugateManually (TestAdjugate):

    @unittest.skipIf(SKIP_TEST, 'Ran 1 test in 261.387s')
    def testConditionNumber(self):
        REPEATS = 8
        with open(f'{OUTDIR}/Python/Output/MatrixCondition.txt', 'w') as fc:
            fc.write("Size\tType\tNoise\tCondition Number"
                    "\tDeterminant Value\tDeterminant Uncertainty\tDeterminant Precision"
                    "\n")
            def write(ssMat, matrixType, noise):
                cond = numpy.linalg.cond(numpy.asarray([[var.value() for var in row] for row in ssMat]))
                det, ssAdj = adjugate(ssMat)
                detUnc = det.uncertainty()
                fc.write(f'{size}\t{matrixType}\t{noise}\t{cond}'
                         f'\t{det.value()}\t{detUnc}\t{detUnc/abs(det.value())}\n')  
                fc.flush()

            for size in range(2, Adjugate.MAX_SIZE):
                ssHilbert = createHilbertMatrix(size)
                for noise in (0, 1e-20, 1e-18, 1e-16, 1e-14):
                    write(addNoise(ssHilbert, noise), 'Hilbert', noise)
                for repeat in range(REPEATS):
                    ssMat = tuple([tuple([VarDbl(random.normalvariate()) for col in range(size)]) 
                                for row in range(size)])
                    write(ssMat, 'Random', 0)
 
    @unittest.skipIf(SKIP_TEST, 'Too slow')
    def testAdjugate(self):
        '''
        '''
        Adjugate.dump()


class TestFFTDumpFile (unittest.TestCase):

    @unittest.skipIf(SKIP_TEST, "After 39 seconds, Found none betwwen [0.1990569896377899, 0.19905698963778992]")
    def test_NotStableException(self):
        sTayler = [-i - 1 if (i & 1) else i + 1 for i in range(Taylor.momentum.maxOrder)]
        lower = 0.199
        upper = 0.200
        Taylor.taylor1d(VarDbl(1, lower), 'lower', sTayler, True, True)
        with self.assertRaises(NotMonotonicException):
            Taylor.taylor1d(VarDbl(1, upper), 'lower', sTayler, True, True)
        for n in range(1000):
            try:
                mid = (lower + upper) * 0.5
                Taylor.taylor1d(VarDbl(1, mid), 'lower', sTayler, True, True)
                lower = mid
                upper = mid
            except NotStableException:
                dumpPath = f'{OUTDIR}/Python/Output/Pow_1_{mid}_-2.txt'      
                with self.assertRaises(NotStableException):
                    Taylor.taylor1d(VarDbl(1, upper), 'lower', sTayler, True, True, dumpPath=dumpPath)
                sInput, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
                self.assertAlmostEqual(sInput['input'].value(), 1)
                self.assertAlmostEqual(sInput['input'].uncertainty(), mid)
                self.assertEqual(len(sExpansion), 223)
                self.assertEqual(sExpansion[-1].monotonics, 25)
                self.assertEqual(out, 'NotStableException')
                break
        else:
            self.fail(f"Found none betwwen [{lower}, {upper}]")


class TestRegressiveSin (unittest.TestCase):

    @unittest.skipIf(SKIP_TEST, "Ran 1 test in 702s")
    def test_18(self):
        sin = RegressiveSin(18)
        sin.calc()



class TestConvergence (unittest.TestCase):
    HEADER = 'X\tEdge\tBias\tValue\tUncertainty\tException\n'

    @staticmethod
    def _search(dxMin: float, dxMax: float, dxRes:int, sX: list[float], 
                varFunc:Callable[[float, float], VarDbl], valFunc:Callable[[float], float],
                path: str):
        with open(path, 'w') as f:
            f.write(TestConvergence.HEADER)
            for x in sX:
                excpt = None
                res = None
                iMin = int(dxMin * dxRes)
                iMax = int(dxMax * dxRes)
                while iMin + 1 < iMax:
                    iMid = int((iMin + iMax)/2)
                    try:
                        dx = iMid / dxRes
                        res = varFunc(x, dx)
                        iMin = iMid
                    except BaseException as ex:
                        excpt = ex
                        iMax = iMid
                if res is None:
                    raise ValueError(f'No result for x={x}, dxMin={dxMin}, dxMax={dxMax}, exception={excpt}')
                if excpt is None:
                    raise ValueError(f'Noexception for x={x}, dxMin={dxMin}, dxMax={dxMax}')
                f.write(f'{x}\t{iMin/dxRes}\t{res.value() - valFunc(x)}\t{res.value()}\t{res.uncertainty()}\t{excpt}\n')
                f.flush()

    @unittest.skipIf(SKIP_TEST, 'Ran 1 test in 2088s')
    def test_pow(self):
        Taylor.pow(VarDbl(1, 0.19929), -1.75, dumpPath=f'{OUTDIR}/Python/Output/Pow_1_0.19929_-1.75.txt')

        DIVIDS = 20
        TestConvergence._search(0.19, 0.21, 100000,
                [i/20 for i in range(-3*DIVIDS, 4*DIVIDS + 1) if (i < 0) or ((i % DIVIDS) != 0)],
                lambda x, dx: Taylor.pow(VarDbl(1, dx), x), lambda x: 1,
                f'{OUTDIR}/Python/Output/PowEdge.txt')

    @unittest.skipIf(SKIP_TEST, 'Ran 1 test in 4187s')
    def test_pow_uniform(self):
        Taylor.pow(VarDbl(1, 0.57), -3, momentum=momentum.UNIFORM,
                   dumpPath=f'{OUTDIR}/Python/Output/Pow_1_0.57_-3.Uniform.txt')
        Taylor.pow(VarDbl(1, 0.58), 2.9, momentum=momentum.UNIFORM,
                   dumpPath=f'{OUTDIR}/Python/Output/Pow_1_0.58_2.9.Uniform.txt')

        DIVIDS = 20
        TestConvergence._search(0.57, 0.59, 100000,
                [i/20 for i in range(-3*DIVIDS, 4*DIVIDS + 1) if (i < 0) or ((i % DIVIDS) != 0)],
                lambda x, dx: Taylor.pow(VarDbl(1, dx), x, momentum=momentum.UNIFORM), lambda x: 1,
                f'{OUTDIR}/Python/Output/PowEdge.Uniform.txt')

    def test_sin(self):
        DIVIDS = 64
        TestConvergence._search(0.3, 5.0, 1000,
                [math.pi*i/DIVIDS for i in range(-1*DIVIDS, 1*DIVIDS + 1)],
                lambda x, dx: Taylor.sin(VarDbl(x, dx)), 
                lambda x: math.sin(x),
                f'{OUTDIR}/Python/Output/SinEdge.txt')
        
    def test_sin_uniform(self):
        DIVIDS = 64
        TestConvergence._search(0.3, 5.0, 1000,
                [math.pi*i/DIVIDS for i in range(-1*DIVIDS, 1*DIVIDS + 1)],
                lambda x, dx: Taylor.sin(VarDbl(x, dx), momentum=momentum.UNIFORM), 
                lambda x: math.sin(x),
                f'{OUTDIR}/Python/Output/SinEdge.Uniform.txt')
        

    def test_exp(self):
        Taylor.exp(VarDbl(0, 19), dumpPath=f'{OUTDIR}/Python/Output/Exp_0_16.txt')
        with self.assertRaises(NotMonotonicException):
            Taylor.exp(VarDbl(0, 20), dumpPath=f'{OUTDIR}/Python/Output/Exp_0_17.txt')

        dumpPath = f'{OUTDIR}/Python/Output/ExpEdge.txt'
        TestConvergence._search(19, 20, 1000,
                (0, 1, -1, 2, -2, 5, -5, 10, -10, 20, -20, 50, -50, 100, -100),
                lambda x, dx: Taylor.exp(VarDbl(x, dx)), 
                lambda x: math.exp(x),
                dumpPath)
        with open(dumpPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestConvergence.HEADER)  
            for line in f:
                sWords = line.split('\t')
                x, edge, bias, val, unc = map(float, sWords[:-1])
                exception = sWords[-1].strip()
                self.assertEqual(edge, 19.864)
                prec = unc / val
                self.assertAlmostEqual(prec, 1681.7672471)
                self.assertTrue(exception.startswith('NotMonotonicException'))


    def test_log(self):
        dumpPath = f'{OUTDIR}/Python/Output/LogEdge.txt'
        TestConvergence._search(0.20, 0.21, 100000,
                (1, 2, 0.5, 5, 0.2, 10, 0.1, 20, 0.05, 50, 0.02, 100, 0.01, 200, 0.005, 500, 0.002, 1000, 0.001),
                lambda x, dx: Taylor.log(VarDbl(x, x*dx)), 
                lambda x: math.log(x),
                dumpPath)
        with open(dumpPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestConvergence.HEADER)  
            for line in f:
                sWords = line.split('\t')
                x, edge, bias, val, unc = map(float, sWords[:-1])
                exception = sWords[-1].strip()
                self.assertEqual(edge, 0.20086)
                self.assertAlmostEqual(unc, 0.2130506)
                self.assertTrue(exception.startswith('NotMonotonicException'))
        
    def test_log_uniform(self):
        dumpPath = f'{OUTDIR}/Python/Output/LogEdge.uniform.txt'
        TestConvergence._search(0.57, 0.59, 100000,
                (1, 2, 0.5, 5, 0.2, 10, 0.1, 20, 0.05, 50, 0.02, 100, 0.01, 200, 0.005, 500, 0.002, 1000, 0.001),
                lambda x, dx: Taylor.log(VarDbl(x, x*dx), momentum=momentum.UNIFORM), 
                lambda x: math.log(x),
                dumpPath)
        with open(dumpPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestConvergence.HEADER)  
            for line in f:
                sWords = line.split('\t')
                x, edge, bias, val, unc = map(float, sWords[:-1])
                exception = sWords[-1].strip()
                self.assertEqual(edge, 0.57899)
                self.assertAlmostEqual(unc, 1.0350349)
                self.assertTrue(exception.startswith('NotMonotonicException'))
        
