'''
Unit test that needs more than one minutes to finish.

To run a particular test:
    VarianceArithmetic\Python> python -m unittest testManual.<class>.<method>
'''
from collections.abc import Callable
import datetime
import logging
import math
import numpy
import os
import random
import os
import scipy
import time
import unittest

from fft import NoiseType, FFT_Order, FFT_Step
from indexSin import SinSource
from matrix import createHilbertMatrix, addNoise
from matrix import adjugate
import momentum
from regressiveSin import RegressiveSin
from taylor import Taylor, NotMonotonicException, NotStableException
from varDbl import VarDbl

from testMatrix import Adjugate, TestAdjugate


SKIP_TEST = False

OUTDIR = f'./Python/Output' if os.getcwd().endswith('VarianceArithmetic') else "./Output"

logger = logging.getLogger(__name__)


MAX_COUNT = 100000

def calcStat(sSample:numpy.array, samples:int) -> tuple[tuple[float, float]]:
    '''
    Reuse the random sampling array {sSample} to calculate the sample mean and standard deviation
    '''
    return tuple([(sSample[i: i+samples].mean(), sSample[i: i+samples].std())
                  for i in range(0, min(len(sSample), MAX_COUNT*samples), samples)])


class TestNormal (unittest.TestCase):
    SAMPLES = (2,3,5,7,10,20,50,100,200,500,1000,2000,5000,10000)
    KAPPAS = (1,1.25,1.5,1.75, 2,2.25,2.5,2.75, 3,3.5, 4, 5, 6)
    MIN_COUNT = 10000

    leakPath = f"{OUTDIR}/NormalLeakage.txt"
    leakHeader = 'Samples\tKappa\tCount\tMean\tDeviation\tRange\n'
    samplesPath = f"{OUTDIR}/NormalSamples.txt"
    samplesHeader = 'Samples\tSlope\tSlopeErr\tIntercept\tInterceptError\tCorrelation\tSamples\n'

    @staticmethod
    def normLeak(mu:float, sigma:float, kapp:float):
        '''
        The sample bounding leakage of normal distribution for the bounding factor {kappa}.
        {mu} is the sample mean, {sigma} is the sample standard deviation.
        '''
        return 1 - scipy.special.erf(abs(mu + kapp*sigma)/math.sqrt(2))*0.5 \
                 - scipy.special.erf(abs(mu - kapp*sigma)/math.sqrt(2))*0.5

    @staticmethod
    def calcNorm(leakPath:str):
        print(f'Start generating {TestNormal.MIN_COUNT} * {TestNormal.SAMPLES[-1]} samples for Normal distribution at {datetime.datetime.now()}')
        sSample = scipy.stats.norm.rvs(size = TestNormal.MIN_COUNT * TestNormal.SAMPLES[-1])
        with open(leakPath, 'w') as f:
            f.write(TestNormal.leakHeader)
            for k in TestNormal.KAPPAS:
                leak = 1 - scipy.special.erf(k/math.sqrt(2))
                factor = scipy.special.erfinv(1 - leak) * math.sqrt(2)
                f.write(f'0\t{k}\t0\t{leak}\t0\t{factor}\n')
            for samples in TestNormal.SAMPLES:
                print(f'Start calculate samples={samples} for Normal distribution at {datetime.datetime.now()}')
                sStat = calcStat(sSample, samples)
                for k in TestNormal.KAPPAS:
                    sLeak = [TestNormal.normLeak(mu, sigma, k) for mu, sigma in sStat]
                    mean = numpy.mean(sLeak)
                    factor = scipy.special.erfinv(1 - mean) * math.sqrt(2)
                    f.write(f'{samples}\t{k}\t{len(sLeak)}\t{mean}\t{numpy.std(sLeak)}\t{factor}\n')
                    f.flush()
        print(f'Finsh writing to {leakPath}')

    def test_normLeak(self):
        for k in range(1, 7):
            self.assertAlmostEqual(TestNormal.normLeak(0, 1, k), 
                    1 - scipy.special.erf(k/math.sqrt(2))/2*2)
            self.assertAlmostEqual(TestNormal.normLeak(0.1, 1, k), 
                    1 - scipy.special.erf((k + 0.1)/math.sqrt(2))/2 
                      - scipy.special.erf((k - 0.1)/math.sqrt(2))/2)
            self.assertAlmostEqual(TestNormal.normLeak(0, 0.9, k), 
                    1 - scipy.special.erf((k*0.9)/math.sqrt(2))/2*2)
            self.assertAlmostEqual(TestNormal.normLeak(0, 1.1, k), 
                    1 - scipy.special.erf((k*1.1)/math.sqrt(2))/2*2)
            self.assertAlmostEqual(TestNormal.normLeak(0.1, 0.9, k), 
                    1 - scipy.special.erf((k*0.9 + 0.1)/math.sqrt(2))/2 
                      - scipy.special.erf((k*0.9 - 0.1)/math.sqrt(2))/2)
            
    @unittest.skipIf(SKIP_TEST, 'Ran 1 tests in 121s')
    def test_calc(self):
        TestNormal.calcNorm(TestNormal.leakPath)
        self.assertTrue(os.path.isfile(TestNormal.leakPath))
        time.sleep(1)
        with open(TestNormal.leakPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestNormal.leakHeader)
            for k in TestNormal.KAPPAS:
                line = next(f)
                n, kappa, cnt, leak, std, fact = map(float, line.split('\t'))
                self.assertEqual(n, 0)
                self.assertEqual(kappa, k)
                self.assertEqual(cnt, 0)
                self.assertAlmostEqual(leak, 1 - scipy.special.erf(k/math.sqrt(2))/2*2)
                self.assertEqual(std, 0)
                self.assertAlmostEqual(k, fact, places=6)
            sLeakPrev = [0] * len(TestNormal.KAPPAS)
            sFactPrev = [0] * len(TestNormal.KAPPAS)
            for samples in TestNormal.SAMPLES:
                for i, k in enumerate(TestNormal.KAPPAS):
                    line = next(f)
                    n, kappa, cnt, leak, std, fact = map(float, line.split('\t'))
                    self.assertEqual(n, samples)
                    self.assertAlmostEqual(kappa, k)
                    self.assertLessEqual(cnt, MAX_COUNT)
                    if i:
                        self.assertGreater(sLeakPrev[i - 1], leak)
                        self.assertLess(sFactPrev[i - 1], fact)
                    if 3 < samples:
                        self.assertGreater(sLeakPrev[i], leak)
                        self.assertLess(sFactPrev[i], fact)
                    sLeakPrev[i] = leak
                    sFactPrev[i] = leak

    def test_samples(self):
        '''
        extract the sample count from 1/kappa.
        '''
        with open(TestNormal.leakPath) as f, open(TestNormal.samplesPath, 'w') as fw:
            fw.write(TestNormal.samplesHeader)
            hdr = next(f)
            self.assertEqual(hdr, TestNormal.leakHeader)
            sX = [1/k for k in TestNormal.KAPPAS]
            sY = []
            for k in TestNormal.KAPPAS:
                line = next(f)
                sY.append(1/float(line.split('\t')[-1]))
            fit = scipy.stats.linregress(sX, sY)
            fw.write(f'0\t{fit.slope}\t{fit.stderr}\t{fit.intercept}\t{fit.intercept_stderr}\t{fit.rvalue}\t{1/fit.intercept}\n')
            for samples in TestNormal.SAMPLES:
                sY = []
                for k in enumerate(TestNormal.KAPPAS):
                    line = next(f)
                    sY.append(1/float(line.split('\t')[-1]))
                fit = scipy.stats.linregress(sX, sY)
                fw.write(f'{samples}\t{fit.slope}\t{fit.stderr}\t{fit.intercept}\t{fit.intercept_stderr}\t{fit.rvalue}\t{1/fit.intercept}\n')



class TestUniform (unittest.TestCase):
    SAMPLES = (2,3,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
    MIN_COUNT = 1000
    leakPath = f'{OUTDIR}/UniformLeakage.txt'
    leakHeader = 'Samples\tCount\tMean\tDeviation\n'

    @staticmethod
    def uniformLeak(mu:float, sigma:float):
        '''
        The sample bounding leakage of uniform distribution between [0,1].
        {mu} is the sample mean, {sigma} is the sample standard deviation.
        '''
        sigma *= math.sqrt(3)
        err = 0
        if mu > sigma:
            err += mu - sigma
        if mu + sigma < 1:
            err += 1 - (mu + sigma)
        return err

    @staticmethod
    def calcUniform(leakHeader:str):
        print(f'Start generating {TestUniform.MIN_COUNT} * {TestUniform.SAMPLES[-1]} samples for uniform distribution at {datetime.datetime.now()}')
        sSample = scipy.stats.uniform.rvs(size = TestUniform.MIN_COUNT * TestUniform.SAMPLES[-1])

        print(f'Start output to {leakHeader}')
        with open(leakHeader, 'w') as f:
            f.write(TestUniform.leakHeader)
            for samples in TestUniform.SAMPLES:
                print(f'Start calculate samples={samples} for uniform distribution at {datetime.datetime.now()}')
                sStat = calcStat(sSample, samples)
                sLeak = [TestUniform.uniformLeak(mu, sigma) for mu, sigma in sStat]
                f.write(f'{samples}\t{len(sLeak)}\t{numpy.mean(sLeak)}\t{numpy.std(sLeak)}\n')
                f.flush()


    def test_uniformLeak(self):
        sigma = 0.5/math.sqrt(3)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.5, sigma), 0)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.4, sigma), 0.1)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.6, sigma), 0.1)
        sigma = 0.6/math.sqrt(3)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.5, sigma), 0)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.6, sigma), 0)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.7, sigma), 0.1)
        sigma = 0.4/math.sqrt(3)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.5, sigma), 0.2)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.6, sigma), 0.2)
        self.assertAlmostEqual(TestUniform.uniformLeak(0.7, sigma), 0.3)

    @unittest.skipIf(SKIP_TEST, 'Ran 1 tests in 460s')
    def test_calc(self):
        TestUniform.calcUniform(TestUniform.leakPath)
        self.assertTrue(os.path.isfile(TestUniform.leakPath))
        time.sleep(1)
        with open(TestUniform.leakPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestUniform.leakHeader)
            for samples in TestUniform.SAMPLES:
                line = next(f)
                n, cnt, leak, std = map(float, line.split('\t'))
                self.assertEqual(n, samples)
                self.assertLessEqual(cnt, MAX_COUNT)
                if 5 <= n:
                    self.assertGreater(leakPrev, leak)
                    self.assertGreater(stdPrev, std)
                leakPrev = leak
                stdPrev = std
 



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
        with open(f'{OUTDIR}/MatrixCondition.txt', 'w') as fc:
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

            for size in range(4, Adjugate.MAX_SIZE):
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
                dumpPath = f'{OUTDIR}/Pow_1_{mid}_-2.txt'      
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
        Taylor.pow(VarDbl(1, 0.19929), -1.75, dumpPath=f'{OUTDIR}/Pow_1_0.19929_-1.75.txt')

        DIVIDS = 20
        TestConvergence._search(0.19, 0.21, 100000,
                [i/20 for i in range(-3*DIVIDS, 4*DIVIDS + 1) if (i < 0) or ((i % DIVIDS) != 0)],
                lambda x, dx: Taylor.pow(VarDbl(1, dx), x), lambda x: 1,
                f'{OUTDIR}/PowEdge.txt')

    @unittest.skipIf(SKIP_TEST, 'Ran 1 test in 4187s')
    def test_pow_uniform(self):
        Taylor.pow(VarDbl(1, 0.57), -3, momentum=momentum.UNIFORM,
                   dumpPath=f'{OUTDIR}/Pow_1_0.57_-3.Uniform.txt')
        Taylor.pow(VarDbl(1, 0.58), 2.9, momentum=momentum.UNIFORM,
                   dumpPath=f'{OUTDIR}/Pow_1_0.58_2.9.Uniform.txt')

        DIVIDS = 20
        TestConvergence._search(0.57, 0.59, 100000,
                [i/20 for i in range(-3*DIVIDS, 4*DIVIDS + 1) if (i < 0) or ((i % DIVIDS) != 0)],
                lambda x, dx: Taylor.pow(VarDbl(1, dx), x, momentum=momentum.UNIFORM), lambda x: 1,
                f'{OUTDIR}/PowEdge.Uniform.txt')

    def test_sin(self):
        DIVIDS = 64
        TestConvergence._search(0.3, 5.0, 1000,
                [math.pi*i/DIVIDS for i in range(-1*DIVIDS, 1*DIVIDS + 1)],
                lambda x, dx: Taylor.sin(VarDbl(x, dx)), 
                lambda x: math.sin(x),
                f'{OUTDIR}/SinEdge.txt')
        
    def test_sin_uniform(self):
        DIVIDS = 64
        TestConvergence._search(0.3, 5.0, 1000,
                [math.pi*i/DIVIDS for i in range(-1*DIVIDS, 1*DIVIDS + 1)],
                lambda x, dx: Taylor.sin(VarDbl(x, dx), momentum=momentum.UNIFORM), 
                lambda x: math.sin(x),
                f'{OUTDIR}/SinEdge.Uniform.txt')
        

    def test_exp(self):
        Taylor.exp(VarDbl(0, 19), dumpPath=f'{OUTDIR}/Exp_0_16.txt')
        with self.assertRaises(NotMonotonicException):
            Taylor.exp(VarDbl(0, 20), dumpPath=f'{OUTDIR}/Exp_0_17.txt')

        dumpPath = f'{OUTDIR}/ExpEdge.txt'
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
        dumpPath = f'{OUTDIR}/LogEdge.txt'
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
        dumpPath = f'{OUTDIR}/LogEdge.uniform.txt'
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
        
