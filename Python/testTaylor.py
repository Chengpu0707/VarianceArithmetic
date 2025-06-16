from collections.abc import Callable
import math
import logging
import numpy as np
import os
import random
import unittest
import sys

from histo import Stat, Histo
from taylor import Taylor, Taylor1dException, NotFiniteException, NotPositiveException, NotMonotonicException, NotStableException
from varDbl import VarDbl, InitException

logger = logging.getLogger(__name__)

taylor = Taylor.default()

class TestConst (unittest.TestCase):

    def test_maxOrder(self):  
        self.assertEqual(taylor.momentum.maxOrder, 448)

    def test_bounding(self):  
        self.assertEqual(taylor.momentum.bounding, 5)

    def test_leakage(self):
        self.assertAlmostEqual(taylor.leakage, 5.7330314e-7)


def _dump_test(self, test:str, func, npfun, sX:tuple[float],
               sDev=[0.2] + [math.pow(10,-err) for err in range(1, 17)] + [0],
               ignoreAssertError:bool=True ):
    HIST_RANGE = 3
    HIST_DIVIDES = 5
    HIST_BINS = 3 * 2 * HIST_DIVIDES
    SAMPLES = 10000
    LOG_PATH = f'./Python/Output/{test}Var.log'
    if os.path.isfile(LOG_PATH):
        os.remove(LOG_PATH)
    logging.basicConfig(filename=LOG_PATH, encoding='utf-8', level=logging.DEBUG,
                        format='%(asctime)s:%(levelname)s:%(message)s')
       
    with open(f'./Python/Output/{test}Var.txt', 'w') as f:
        f.write(f"NoiseType\tNoise\tX\t{test}\tError Deviation\tError Minimum\tError Maximum"
                "\tValue Deviation\tUncertainty\tMean\tBias")
        for i in range(HIST_BINS):
            f.write(f'\t{(i - HIST_BINS//2 + 0.5) /HIST_DIVIDES:.1f}')
        f.write('\n')

        for dev in sDev:
            ssNoise = {'Gaussian': np.random.normal(0, dev, SAMPLES),
                       'Uniform': np.array([(i * 2 - SAMPLES) / SAMPLES * dev * math.sqrt(3) for i in range(SAMPLES + 1)])}
            for x in sX:
                base = func(x)
                try:
                    var = self.validate(x, dev, exception=AssertionError if ignoreAssertError else None)
                except AssertionError as ex:
                    logger.info(f'Ignore {test}({x}+/-{dev}) due to AssertionError {ex}')
                    continue
                except Taylor1dException as ex:
                    logger.info(f'Ignore {test}({x}+/-{dev}) due to {ex}')
                    continue
                except BaseException as ex:
                    logger.info(f'Ignore {test}({x}+/-{dev}) due to BaseException {ex}')
                    continue
                if var.uncertainty() <= 0:
                    continue
                
                for NoiseType, sNoise in ssNoise.items():
                    try:
                        sError = npfun(x, sNoise) - base
                    except Exception as ex:
                        logger.info(f'Ignore {test}({x}+/-{dev}) due to numpy Exception {ex}')
                        continue
                    except BaseException:
                        logger.info(f'Ignore {test}({x}+/-{dev}) due to numpy BaseException {ex}')
                        continue
                    if not math.isfinite(np.min(sError)):
                        logger.info(f'Ignore {test}({x}+/-{dev}) due to numpy nan')
                        continue
                    valDev = np.std(sError)
                    sNorm = sError / var.uncertainty()
                    f.write(f"{NoiseType}\t{dev}\t{x}\t{base}\t{np.std(sNorm)}\t{np.min(sNorm)}\t{np.max(sNorm)}")
                    f.write(f"\t{valDev}\t{var.uncertainty()}\t{np.mean(sError)}\t{var.value() - base}")
                    hist, bin_edges = np.histogram(sNorm, bins=HIST_BINS, range=(-HIST_RANGE, +HIST_RANGE), density=True)
                    for i in range(HIST_BINS):
                        f.write(f'\t{hist[i]}')
                    f.write('\n')

def _validate(self, func, arg, uncertainty, exception, 
              dumpPath, lines = 0):
    if dumpPath:
        try:
            os.remove(dumpPath)
        except:
            pass
    try:
        res = func(arg, uncertainty, dumpPath)
        if exception is not None:
            self.fail(f'No {exception} for {func.__name__}({arg}, {uncertainty})={res} not throw')
        if dumpPath:
            self.assertTrue(os.path.isfile(dumpPath))
            with open(dumpPath) as f:
                cnt = 0
                for line in f:
                    cnt += 1
                if lines:
                    self.assertEqual(lines, cnt)
                sVal = list(map(float, line.split('\t')))
                self.assertAlmostEqual(res.value(), sVal[0])
                self.assertAlmostEqual(res.variance(), sVal[1] + sVal[2])
        return res
    except BaseException as ex:
        if (exception is None) or (not isinstance(ex, exception)): 
            raise ex
        if dumpPath:
            self.assertTrue(os.path.isfile(dumpPath))
            with open(dumpPath) as f:
                cnt = 0
                for line in f:
                    cnt += 1
                if lines:
                    self.assertEqual(lines, cnt)
                self.assertEqual(f'{exception.__name__}', line.strip().split("\t")[0])

    

class TestExp (unittest.TestCase):
    @staticmethod
    def exp(arg, uncertainty, dumpPath):
        var = VarDbl(arg, uncertainty)
        return taylor.exp(var, dumpPath)

    def validate(self, exp, uncertainty, exception=None, 
                 valueDelta=2e-5, varianceDelta=2e-5, dumpPath=None, lines=0) -> VarDbl:
        res = _validate(self, TestExp.exp, exp, uncertainty, 
                        exception = None if (exception == InitException) or (exception == AssertionError) else exception, 
                        dumpPath=dumpPath, lines=lines)
        if res is None:
            return
        try:
            prec = res * math.exp(-exp) - 1
        except InitException as ex:
            if (exception is None) or (not isinstance(ex, InitException)): 
                raise ex
            else:
                return res
        try:
            delta = math.pow(uncertainty, 2)/2 + math.pow(uncertainty, 4)/8 +\
                    math.pow(uncertainty, 6)/48 + math.pow(uncertainty, 8)/384
            self.assertAlmostEqual(prec.value()/delta, 1, delta=valueDelta)
            unc = math.pow(uncertainty, 2) + math.pow(uncertainty, 4)*3/2 +\
                  math.pow(uncertainty, 6)*7/6 + math.pow(uncertainty, 8)*5/8
            self.assertAlmostEqual(prec.variance() / unc, 1, delta=varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res

    def test_range(self):
        self.assertAlmostEqual(709.78, math.log(sys.float_info.max), delta=0.1)

        self.validate(0, 0.1, 
                dumpPath='./Python/Output/exp_1_0.1.txt', lines=168)
        self.validate(0, 0.2, varianceDelta=3e-5)
        self.validate(0, 1, valueDelta=4e-4, varianceDelta=9e-2)

    def test_limit_1_10th(self):
        self.validate(0.1, 0.1)
        self.validate(0.01, 0.1)

        self.validate(-391, 0.1, exception=InitException)     # prec has inf variance
        self.validate(-390, 0.1, varianceDelta=1)   # loss of variance due to small value  
        self.validate(-200, 0.1)

        self.validate(195, 0.1, exception=InitException)      # res has inf variance
        self.validate(194, 0.1)
        self.validate(100, 0.1)
 
    def test_limit_2_10th(self):
        self.validate(0, 0.2, varianceDelta=3e-5)

        self.validate(-392, 0.2, exception=InitException)     # prec has inf variance
        self.validate(-391, 0.2, varianceDelta=1)   # loss of variance due to small value
        self.validate(-200, 0.2, varianceDelta=3e-5)

        self.validate(196, 0.2, exception=NotFiniteException)
        self.validate(195, 0.2, valueDelta=3e-2, varianceDelta=4e-2)
        self.validate(100, 0.2, varianceDelta=3e-5)
 
    def test_limit_1(self):
        self.validate(0, 1, valueDelta=4e-4, varianceDelta=9e-2)

        self.validate(-392, 1, exception=InitException)       # prec has inf variance
        self.validate(-391, 1, valueDelta=4e-4, varianceDelta=1)
        self.validate(-200, 1, valueDelta=4e-4, varianceDelta=9e-2)

        self.validate(196, 1, exception=NotFiniteException)
        self.validate(195, 1, valueDelta=7e-1, varianceDelta=3)
        self.validate(100, 1, valueDelta=4e-4, varianceDelta=9e-2)
 
    def test_exp_min_input_precision(self):
        for x in (0, 1, -1, 10, -10, 100, -100):
            RES = 10000
            iMax = RES
            iMin = RES // 2
            while iMin + 1 < iMax:
                iMid = int((iMin + iMax)/2)
                res = taylor.exp(VarDbl(x, iMid/RES))
                if res.value() < res.uncertainty():
                    iMax = iMid
                else:
                    iMin = iMid
            self.assertEqual(iMid, 8328)


    @staticmethod
    def func(x):
        return math.exp(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.exp(sNoise + x)
 
    def test_dump(self):
        _dump_test(self, 'exp', TestExp.func, TestExp.npfunc,
                  (-100, -50, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 50, 100),
                  sDev = [math.pow(10,-err) for err in range(0, 17)])


class TestLog (unittest.TestCase):
    @staticmethod
    def log(arg, uncertainty, dumpPath):
        var = VarDbl(arg, uncertainty)
        return taylor.log(var, dumpPath)

    def validate(self, x, uncertainty, 
                 exception=None, valueDelta=2e-2, varianceDelta=2e-2,
                 dumpPath=None, lines=0) -> VarDbl:
        res = _validate(self, TestLog.log, x, uncertainty, 
                        exception = None if (exception == AssertionError) else exception, 
                        dumpPath=dumpPath, lines=lines)
        if res is None:
            return
        try:
            precOut = res - math.log(x)
            precIn = abs(uncertainty/x)
            bias = - math.pow(precIn, 2)/2 - math.pow(precIn, 4)/4 - math.pow(precIn, 6)/6 - math.pow(precIn, 8)/8
            self.assertAlmostEqual(precOut.value() / bias, 1, delta=valueDelta)
            unc = math.pow(precIn, 2) + math.pow(precIn, 4)*9/8 + math.pow(precIn, 6)*119/24 + math.pow(precIn, 8)*991/32
            self.assertAlmostEqual(precOut.variance() /unc, 1, delta=varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res

    def test_exception(self):
        self.assertAlmostEqual(709.78, math.log(sys.float_info.max), delta=0.1)

        self.validate(1, 0.1,
                dumpPath='./Python/Output/log_1_0.1.txt', lines=168)
        self.validate(2, 0.2)

        self.validate(0.5, 0.101, exception=NotMonotonicException)
        self.validate(0.5, 0.100, valueDelta=0.05, varianceDelta=0.07)

        self.validate(1, 0.201, exception=NotMonotonicException)
        self.validate(1, 0.200, valueDelta=0.05, varianceDelta=0.07)

    def test_max_out_uncertainty(self):
        for x in (0.5, 1, 2, 5, 10, 20, 50):
            res = taylor.log(VarDbl(x, x/5))
            self.assertAlmostEqual(res.uncertainty(), 0.2120048)

    @staticmethod
    def func(x):
        return math.log(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.log(sNoise + x)
 
    def test_dump(self):
        _dump_test(self, 'log', TestLog.func, TestLog.npfunc,
                  (1/32, 1/20, 1/16, 0.1, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                  ignoreAssertError=True)
        

class TestSin (unittest.TestCase):
    @staticmethod
    def sin(arg, uncertainty, dumpPath):
        var = VarDbl(arg, uncertainty)
        return taylor.sin(var, dumpPath)

    def validate(self, x, uncertainty, exception = None, 
                 valueDelta = 2e-5, varianceDelta = 2e-4, 
                 dumpPath = None, lines = 0) -> VarDbl:
        res = _validate(self, TestSin.sin, x, uncertainty, 
                        exception = None if (exception == AssertionError) else exception, 
                        dumpPath=dumpPath, lines=lines)
        if res is None:
            return
        try:
            prec = res - math.sin(x)
            bias = math.sin(x) * (-math.pow(uncertainty, 2)/2 + math.pow(uncertainty, 4)/8 +\
                  -math.pow(uncertainty, 6)/48 + math.pow(uncertainty, 8)/384)
            if bias:
                self.assertAlmostEqual(prec.value() / bias, 1, delta = valueDelta)
            else:
                self.assertAlmostEqual(prec.value(), bias, delta = valueDelta * prec.value())
            cos2 = math.cos(x)**2
            var = math.pow(uncertainty, 2) * cos2 - math.pow(uncertainty, 4)*(3/2 * cos2 - 1/2) +\
                  math.pow(uncertainty, 6)*(7/6 * cos2 - 1/2) - math.pow(uncertainty, 8)*(5/8 * cos2 - 7/24)
            if var:
                self.assertAlmostEqual(prec.variance() / var, 1, delta = varianceDelta)
            else:
                self.assertAlmostEqual(prec.variance(), var, delta = varianceDelta * prec.variance())
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res
        except BaseException as ex:
            if (exception is None) or (not isinstance(ex, exception)): 
                raise ex
        
    def test_exception_pi(self):
        self.validate(0, 0.1,
                dumpPath='./Python/Output/sin_0_0.1.txt')

        self.validate(0, math.pi/8)
        self.validate(0, math.pi/4, varianceDelta=3e-2)
        self.validate(0, math.pi/2, NotPositiveException,
                dumpPath='./Python/Output/sin_0_0.5.txt')

        self.validate(math.pi, math.pi/8)
        self.validate(math.pi, math.pi/4, valueDelta=8e-5, varianceDelta=3e-2)
        self.validate(math.pi, math.pi/2, NotPositiveException)

        res1 = self.validate(- math.pi/64, 0.1)
        res2 = self.validate(+ math.pi/64, 0.1)
        self.assertAlmostEqual(-res1.value(), res2.value())
        self.assertAlmostEqual(res1.uncertainty(), res2.uncertainty())

    def test_exception_half_pi(self):
        self.validate(math.pi/2, 0.1)
        self.validate(math.pi/2, math.pi/8, varianceDelta=2e-3)
        res = self.validate(math.pi/2, math.pi/4, valueDelta=8e-5, varianceDelta=8e-2)
        self.assertAlmostEqual(res.value() - math.sin(math.pi/2), -0.2653961)
        self.validate(math.pi/2, math.pi/2, NotPositiveException)

        self.validate(-math.pi/2, 0.1)
        self.validate(-math.pi/2, math.pi/8, varianceDelta=2e-3)
        res = self.validate(-math.pi/2, math.pi/4, valueDelta=8e-5, varianceDelta=8e-2)
        self.assertAlmostEqual(res.value() - math.sin(-math.pi/2), 0.2653961)
        self.validate(-math.pi/2, math.pi/2, NotPositiveException)

        res1 = self.validate(math.pi/2 - math.pi/64, 0.1)
        res2 = self.validate(math.pi/2 + math.pi/64, 0.1)
        self.assertAlmostEqual(res1.value(), res2.value())
        self.assertAlmostEqual(res1.uncertainty(), res2.uncertainty())

    def test_exception_quater_pi(self):
        res = self.validate(math.pi/4, math.pi/4, valueDelta=8e-5, varianceDelta=2e-3)
        self.assertAlmostEqual(res.value() - math.sin(math.pi/4), -0.1876634)
        self.validate(math.pi/4, math.pi/2, NotPositiveException)

        self.validate(-math.pi/4, 0.1)
        self.validate(-math.pi/4, math.pi/8)
        res = self.validate(-math.pi/4, math.pi/4, valueDelta=8e-5, varianceDelta=2e-3)
        self.assertAlmostEqual(res.value() - math.sin(-math.pi/4), 0.1876634)
        self.validate(-math.pi/4, math.pi/2, NotPositiveException)

        res1 = self.validate(math.pi/4 - math.pi/64, 0.1)
        res2 = self.validate(math.pi/4 + math.pi/64, 0.1)
        self.assertLess(res1.value(), res2.value())
        self.assertGreater(res1.uncertainty(), res2.uncertainty())


    @staticmethod
    def func(x):
        return math.sin(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.sin(sNoise + x)
 
    def test_dump(self):
        _dump_test(self, 'sin', TestSin.func, TestSin.npfunc,
                  np.array([i/16 for i in range(-16, 17)]) * math.pi,
                  sDev = [math.pow(10,-err) for err in range(1, 17)] + [0, 0.2, 0.5])



class TestPow (unittest.TestCase):
    @staticmethod
    def pow(arg, uncertainty, dumpPath):
        var = VarDbl(1, uncertainty)
        return taylor.pow(var, arg, dumpPath)

    def validate(self, exp, uncertainty, exception=None, 
                 valueDelta=2e-5, varianceDelta=2e-5, 
                 dumpPath=None, lines=0) -> VarDbl:
        try:
            res = _validate(self, TestPow.pow, exp, uncertainty, 
                    exception = None if (exception == InitException) or (exception == AssertionError) else exception, 
                    dumpPath=dumpPath, lines=lines)
            if res is None:
                return
            prec = res - 1
            bias = math.pow(uncertainty, 2) * exp*(exp-1)/2 + \
                   math.pow(uncertainty, 4) * exp*(exp-1)*(exp-2)*(exp-3)/24 +\
                   math.pow(uncertainty, 6) * exp*(exp-1)*(exp-2)*(exp-3)*(exp-4)*(exp-5)/720
            if bias:
                self.assertAlmostEqual(prec.value() / bias, 1, delta = valueDelta)
            else:
                self.assertAlmostEqual(prec.value(), bias, delta = prec.value() * valueDelta)
            var = math.pow(uncertainty, 2) * exp*exp + \
                  math.pow(uncertainty, 4) * exp*exp*(exp-1)*(exp-5/3)*3/2 +\
                  math.pow(uncertainty, 6) * exp*exp*(exp-1)*(exp-2)*(exp-2)*(exp-16/7)*7/6
            if var:
                self.assertAlmostEqual(prec.variance() / var, 1, delta = varianceDelta)
            else:
                self.assertAlmostEqual(prec.variance(), var, delta = prec.variance() * varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (exception != AssertionError): 
                raise ex
            else:
                return res
        except BaseException as ex:
            if (exception is None) or (not isinstance(ex, exception)): 
                raise ex
            
    def test_int_exp(self):
        self.validate(0, 0) 
        self.validate(0, 1) 
        self.validate(1, 0) 
        self.validate(1, 1) 
        self.validate(2, 0) 
        self.validate(2, 1, varianceDelta=8e-5) 

    def test_zero(self):
        with self.assertRaises(ValueError) as ex:
            taylor.pow(VarDbl(0, 0.1), -1,
                    dumpPath='./Python/Output/pow_0_0.1_-1.txt')
        self.assertEqual(str(ex.exception), 'math domain error')   # 0^{-1}

        with self.assertRaises(NotFiniteException):
            taylor.pow(VarDbl(0.1, 0.1), -1,
                    dumpPath='./Python/Output/pow_0.1_0.1_-1.txt')

    def test_two(self):
        self.validate(2, 0.1,
                dumpPath='./Python/Output/pow_1_0.1_2.txt')
        self.validate(2 - 1e-6, 0.1, valueDelta=4e-3)  
        self.validate(2 + 1e-6, 0.1, valueDelta=4e-3)  
        self.validate(2, 1, varianceDelta=8e-5,     # due to momentum(4) != 3
                dumpPath='./Python/Output/pow_1_1_2.txt')
        
    def test_inverse(self):
        self.validate(-1, 0.100, valueDelta = 3e-2, varianceDelta = 7e-4,
                      dumpPath="./Python/Output/pow_1_0.1_-1.txt") 
        self.validate(-2, 0.100, valueDelta = 4e-2, varianceDelta = 4e-3,
                      dumpPath="./Python/Output/pow_1_0.1_-2.txt") 
        
        self.validate(-1, 0.201, NotMonotonicException, 
                      dumpPath="./Python/Output/pow_1_0.201_-1.txt") 
        self.validate(-2, 0.200, NotMonotonicException, 
                      dumpPath="./Python/Output/pow_1_0.200_-2.txt") 
        
    def test_continuity(self):
        RES = 10000
        MIN = int(0.2 * RES)
        MAX = int(0.21 * RES)
        with open("./Python/Output/pow_continuity.txt", 'w') as f:
            f.write('n\tdelta n\tUpper Bound\tValue\tUncertainty\tDiff Value\tDiff Uncertainty\tException\n')
            for n in range(10):
                for dx in (1e-6, -1e-6):
                    excpt = None
                    iMin = MIN
                    iMax = MAX
                    while iMin + 1 < iMax:
                        iMid = int((iMin + iMax)/2)
                        try:
                            res = taylor.pow(VarDbl(1, iMid / RES), n + dx)
                            iMin = iMid
                        except BaseException as ex:
                            iMax = iMid
                            excpt = ex
                    upper = iMin/RES
                    std = taylor.pow(VarDbl(1, upper), n)
                    f.write(f'{n}\t{dx}\t{upper}\t{std.value()}\t{std.uncertainty()}\t{res.value() - std.value()}\t{res.uncertainty() - std.uncertainty()}\t{excpt}\n')
                    f.flush()
 

    @staticmethod
    def func(x):
        return 1.0

    @staticmethod
    def npfunc(x, sNoise):
        return np.power(sNoise + 1, x)
 
    def test_dump(self):
        _dump_test(self, 'pow', TestPow.func, TestPow.npfunc,
                np.array([i/10 for i in range(-20, 31, 2)] +[-0.1, 0.1, -0.01, 0.01, -1e-3, 1e-6]),
                sDev = [0.2, 0.197] + [math.pow(10,-err) for err in range(1, 17)])


class TestLibError (unittest.TestCase):

    def testExpLogNoError(self):
        for x in (2, math.sqrt(2), math.pi):
            var = taylor.exp(VarDbl(x))
            res = taylor.log(var) - x
            self.assertEqual(res.value(), 0)

            var = taylor.log(VarDbl(x))
            res = taylor.exp(var) - x
            self.assertEqual(res.value(), 0)
    
    def testExpLog(self):
        with open(f'./Python/Output/ExpLogError.txt', 'w') as f:
            f.write('X\tType\tError\tVarDbl Error\tUncertainty\n')
            for i in range(-200, 201):
                x = i / 100
                var = taylor.log(taylor.exp(VarDbl(x))) - x
                lib = math.log(math.exp(x)) - x
                f.write(f'{x}\tlog(exp(x))\t{lib}\t{var.value()}\t{var.uncertainty()}\n')
                if x <= 0:
                    continue
                var = taylor.exp(taylor.log(VarDbl(x))) - x
                lib = math.exp(math.log(x)) - x
                f.write(f'{x}\texp(log(x))\t{lib}\t{var.value()}\t{var.uncertainty()}\n')




    def testPower(self):
        with open(f'./Python/Output/PowerError.txt', 'w') as f:
            f.write('X\tExp\tVarDbl Error\tVarDbl Uncertainty\tLib Error\tLib Uncertainty\n')
            for i in range(1, 100):
                x = VarDbl(i/10)
                for j in range(-20, 21):
                    if not j:
                        continue
                    try:
                        exp = j / 10
                        var = taylor.pow(taylor.pow(x, exp), 1/exp) - x
                        lib = VarDbl(math.pow(math.pow(x.value(), exp), 1/exp)) - x.value()
                        f.write(f'{x.value()}\t{exp}\t{var.value()}\t{var.uncertainty()}\t{lib.value()}\t{lib.uncertainty()}\n')
                    except BaseException as ex:
                        print(f'x={x} exp ={exp}: {ex}')
                        continue

    def testSin(self):
        dev = 0.01
        with open(f'./Python/Output/SinError_{dev}.txt', 'w') as f:
            f.write('Noise\tLib Error\n')
            x = math.pi/2
            SAMPLES = 100
            for i in range(SAMPLES + 1):
                noise = (i * 2 - SAMPLES) / SAMPLES * dev * math.sqrt(3)
                f.write(f'{noise}\t{math.sin(x + noise) - math.sin(x)}\n')  


class TestDumpFile (unittest.TestCase):

    def test_normal(self):
        dumpPath='./Python/Output/Pow_1_0.197_-2.txt'
        res = taylor.pow(VarDbl(1, 0.197), -2, dumpPath=dumpPath)
        x, sTaylor, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
        self.assertAlmostEqual(res.value(), out.value())
        self.assertAlmostEqual(res.variance(), out.variance())
        self.assertAlmostEqual(x.value(), 1)
        self.assertAlmostEqual(x.uncertainty(), 0.197)
        self.assertEqual(len(sTaylor), taylor.momentum.maxOrder)
        self.assertEqual(len(sExpansion), taylor.momentum.maxOrder//2 - 1)
        self.assertEqual(sExpansion[-1].monotonics, 167)


    def test_NotStableException(self):
        s1dTaylor = [i * (1 if (i & 1) == 0 else -1) for i in range(60)]
        dumpPath='./Python/Output/Pow_1_-2_NotStable.txt'
        with self.assertRaises(NotStableException):
            taylor.taylor1d(VarDbl(1, 0.193), 'test_NotStableException', 
                            s1dTaylor, True, True, maxOrder=60, dumpPath=dumpPath)
        x, sTaylor, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
        self.assertEqual(out, "NotStableException")
        self.assertAlmostEqual(x.value(), 1)
        self.assertAlmostEqual(x.uncertainty(), 0.193)
        self.assertListEqual(sTaylor, s1dTaylor)
        self.assertEqual(len(sExpansion), 29)
        self.assertEqual(sExpansion[-1].monotonics, 28)

    def test_NotMonotonicException(self):
        dumpPath='./Python/Output/Pow_1_0.2_-2.txt'
        with self.assertRaises(NotMonotonicException):
            taylor.pow(VarDbl(1, 0.2), -2, dumpPath=dumpPath)
        x, sTaylor, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
        self.assertEqual(out, "NotMonotonicException")
        self.assertAlmostEqual(x.value(), 1)
        self.assertAlmostEqual(x.uncertainty(), 0.2)
        self.assertEqual(len(sTaylor), taylor.momentum.maxOrder)
        self.assertEqual(len(sExpansion), taylor.momentum.maxOrder//2 - 1)
        self.assertEqual(sExpansion[-1].monotonics, 0)

    def test_NotMonotonicException_AfterMax(self):
        dumpPath = f'./Python/Output/Pow_1_0.19902_2.95.txt'
        res = taylor.pow(VarDbl(1, 0.19902), 2.95, dumpPath=dumpPath)
        x, sTaylor, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
        self.assertAlmostEqual(x.value(), 1)
        self.assertAlmostEqual(x.uncertainty(), 0.19902)
        self.assertEqual(res.value(), out.value())
        self.assertEqual(res.variance(), out.variance())
        self.assertEqual(len(sTaylor), taylor.momentum.maxOrder)
        self.assertEqual(len(sExpansion), taylor.momentum.maxOrder//2 - 1)
        self.assertAlmostEqual(sExpansion[-3].newVar.value(), 9.1629E-20, delta=1e-23)
        self.assertAlmostEqual(sExpansion[-2].newVar.value(), 9.3954E-20, delta=1e-23)
        self.assertAlmostEqual(sExpansion[-1].newVar.value(), 7.6975E-21, delta=1e-24)
        self.assertEqual(sExpansion[-3].monotonics, 220)
        self.assertEqual(sExpansion[-2].monotonics, 220)
        self.assertEqual(sExpansion[-1].monotonics, 221)


    def test_NotPositiveException(self):
        dumpPath = './Python/Output/Sin_0.5_0.3185.txt'
        with self.assertRaises(NotPositiveException):
            taylor.sin(VarDbl(0.5*math.pi, 0.3185*math.pi), dumpPath=dumpPath)
        x, sTaylor, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
        self.assertEqual(out, 'NotPositiveException')
        self.assertAlmostEqual(x.value(), 0.5*math.pi)
        self.assertAlmostEqual(x.uncertainty(), 0.3185*math.pi)
        self.assertEqual(len(sTaylor), taylor.momentum.maxOrder)
        self.assertAlmostEqual(sTaylor[0], 1)
        self.assertAlmostEqual(sTaylor[1], 0)
        self.assertAlmostEqual(sTaylor[2], -0.5)
        self.assertEqual(len(sExpansion), 3)
        self.assertLess(sExpansion[-1].var, 0)

    def test_NotFiniteException(self):
        dumpPath='./Python/Output/pow_1_1_-2.txt'
        with self.assertRaises(NotFiniteException):
            taylor.pow(VarDbl(1, 1), -2, dumpPath=dumpPath)
        x, sTaylor, sExpansion, out = Taylor.verifyDumpFile(self, dumpPath)
        self.assertEqual(out, 'NotFiniteException')
        self.assertAlmostEqual(x.value(), 1)
        self.assertAlmostEqual(x.uncertainty(), 1)
        self.assertEqual(len(sTaylor), taylor.momentum.maxOrder)
        self.assertEqual(len(sExpansion), 122)
        self.assertEqual(sExpansion[-1].monotonics, 0)


class TestConvergence (unittest.TestCase):
    HEADER = 'X\tEdge\tBias\tValue\tUncertainty\tPrecision\tException\n'

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
                f.write(f'{x}\t{iMin/dxRes}\t{res.value() - valFunc(x)}\t{res.value()}\t{res.uncertainty()}\t{res.precision()}\t{excpt}\n')
                f.flush()

    def test_pow(self):
        TestConvergence._search(0.19, 0.21, 100000,
                [i/20 for i in range(-3*20, 4*20 + 1) if (i % 20) != 0],
                lambda x, dx: taylor.pow(VarDbl(1, dx), x), lambda x: 1,
                './Python/Output/PowAtOneEdge.txt')

    def test_sin(self):
        TestConvergence._search(0.3, 0.55, 1000,
                [i/64 for i in range(-1*64, 1*64 + 1)],
                lambda x, dx: taylor.sin(VarDbl(x *math.pi, dx *math.pi)), 
                lambda x: math.sin(x *math.pi),
                './Python/Output/SinEdge.txt')
        

    def test_exp(self):
        dumpPath = './Python/Output/ExpEdge.txt'
        TestConvergence._search(19, 20, 1000,
                (0, 1, -1, 2, -2, 5, -5, 10, -10, 20, -20, 50, -50),
                lambda x, dx: taylor.exp(VarDbl(x, dx)), 
                lambda x: math.exp(x),
                dumpPath)
        with open(dumpPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestConvergence.HEADER)  
            for line in f:
                sWords = line.split('\t')
                x, edge, bias, val, unc, prec = map(float, sWords[:-1])
                exception = sWords[-1].strip()
                self.assertEqual(edge, 19.864)
                self.assertAlmostEqual(prec, 1681.7672471, delta=1e-6)
                self.assertTrue(exception.startswith('NotMonotonicException'))
        
    def test_log(self):
        dumpPath = './Python/Output/LogEdge.txt'
        TestConvergence._search(0.20, 0.21, 100000,
                (1, 2, 0.5, 5, 0.2, 10, 0.1, 20, 0.05, 50, 0.02, 100, 0.01, 200, 0.005, 500, 0.002, 1000, 0.001),
                lambda x, dx: taylor.log(VarDbl(x, x*dx)), 
                lambda x: math.log(x),
                dumpPath)
        with open(dumpPath) as f:
            hdr = next(f)
            self.assertEqual(hdr, TestConvergence.HEADER)  
            for line in f:
                sWords = line.split('\t')
                x, edge, bias, val, unc, prec = map(float, sWords[:-1])
                exception = sWords[-1].strip()
                self.assertEqual(edge, 0.20086)
                self.assertAlmostEqual(unc, 0.2130505830, delta=1e-9)
                self.assertTrue(exception.startswith('NotMonotonicException'))
        


class TestStat (unittest.TestCase):
    @staticmethod
    def writePowerHeader(f, divids:int=5, devs:int=3):
        f.write('Exponent\tInput Value\tInput Uncertainty'
                '\tValue\tUncertainty\tMean\tDeviation'
                '\tNormalized Error Mean\tNormalized Error Deviation\tLess\tMore')
        histo = Histo(divids, devs)
        for bucket in histo.buckets():
            f.write(f'\t{bucket:.1f}')
        f.write('\n')

    @staticmethod
    def calc_power(x:float, dx:float, exp:float, f):
        if dx <= 0:
            return
        CNT = 100000
        try:
            var = VarDbl(x, dx)**exp
        except BaseException as ex:
            raise ex
        varUnc = var.uncertainty()
        stat = Stat()
        histo = Histo(5, 3)
        for j in range(CNT):
            val = random.gauss(x, dx)**exp
            try:
                stat.accum( val )
            except BaseException as ex:
                raise ex
            histo.accum( (val - var.value())/varUnc )
        f.write(f'{exp}\t{x}\t{dx}\t{var.value()}\t{var.uncertainty()}'
                f'\t{stat.mean()}\t{stat.dev()}\t{histo.stat().mean()}\t{histo.stat().dev()}\t{histo.less()}\t{histo.more()}')
        cnt = stat.count() - histo.less() - histo.more()
        for h in histo.histogram():
            f.write(f'\t{h/cnt}')
        f.write('\n')


    def test_square_at_zero(self):
        with open('./Python/Output/SquareAtZero.txt', 'w') as f:
            TestStat.writePowerHeader(f)
            for x in (-0.2, -0.1, 0, 0.1, 0.2):
               TestStat.calc_power(x, 0.2, 2, f)

    def test_natural_number_at_zero(self):
        with open('./Python/Output/NaturalAtZero.txt', 'w') as f:
            TestStat.writePowerHeader(f)
            TestStat.calc_power(0, 0.2, 2, f)
            TestStat.calc_power(0.2, 0.2, 2, f)
            TestStat.calc_power(0, 0.2, 3, f)                
            TestStat.calc_power(0.2, 0.2, 3, f)                
            TestStat.calc_power(-0.2, 0.2, 3, f)                

    def test_square_at_one(self):
        with open('./Python/Output/SquareAtOne.txt', 'w') as f:
            TestStat.writePowerHeader(f)
            for c in (2, 2 - 1e-6, 2 + 1e-6):
               TestStat.calc_power(1, 0.2, c, f)

    def test_invesion_at_one(self):
        EXP = -1
        with self.assertRaises(NotMonotonicException):
            VarDbl(1, 0.2001)**EXP
        res = VarDbl(1, 0.2000)**EXP
        self.assertAlmostEqual(res.value(), 1.0462500, places=6)
        self.assertAlmostEqual(res.uncertainty(), 0.2548432, places=6)

        with open('./Python/Output/InversionAtOne.txt', 'w') as f:
            TestStat.writePowerHeader(f)
            for dx in (0.2, 0.1, 1e-2, 1e-3, 1e-4):
               TestStat.calc_power(1, dx, EXP, f) 

    def test_square_root_at_one(self):
        EXP = 0.5
        with self.assertRaises(NotMonotonicException):
            VarDbl(1, 0.2012)**EXP
        res = VarDbl(1, 0.2011)**EXP
        self.assertAlmostEqual(res.value(), 0.9947278, places=6)
        self.assertAlmostEqual(res.uncertainty(), 0.1025509, places=6)

        with open('./Python/Output/SquareRootAtOne.txt', 'w') as f:
            TestStat.writePowerHeader(f)
            for dx in (0.2011, 1e-1, 1e-2, 1e-3):
               TestStat.calc_power(1, dx, EXP, f) 


class TestRoundingError (unittest.TestCase):
    '''
    small uncertainty of rounding error caues very small expansion order
    '''

    def test_inverse(self):
        inv = taylor.pow(VarDbl(0.1), -1, dumpPath='./Python/Output/Pow_0.1_-1.txt')
        self.assertAlmostEqual(inv.value() / 10, 1)
        self.assertAlmostEqual(inv.uncertainty() / 8.01228266906342e-16, 1)
        self.assertAlmostEqual(math.ulp(inv.value())/1.7763568394002505e-15, 1)
        self.assertAlmostEqual(VarDbl.ulp(inv.value())/1.0255800994045676e-15, 1)

        inv = taylor.pow(VarDbl(0.01), -1, dumpPath='./Python/Output/Pow_0.01_-1.txt')
        self.assertAlmostEqual(inv.value() / 99.99999999999991, 1)
        self.assertAlmostEqual(inv.uncertainty() / 1.0015353336329276e-14, 1)
        self.assertAlmostEqual(math.ulp(inv.value())/1.4210854715202004e-14, 1)
        self.assertAlmostEqual(VarDbl.ulp(inv.value())/8.20464079523654e-15, 1)

        inv = taylor.pow(VarDbl(0.001), -1, dumpPath='./Python/Output/Pow_0.001_-1.txt')
        self.assertAlmostEqual(inv.value() / 1000, 1)
        self.assertAlmostEqual(inv.uncertainty() / 1.2519191670411594e-13, 1)
        self.assertAlmostEqual(math.ulp(inv.value())/1.1368683772161603e-13, 1)
        self.assertAlmostEqual(VarDbl.ulp(inv.value())/6.563712636189232e-14, 1)


class TestImpreciseCoeff (unittest.TestCase):
     
     def test_inverse(self):
        with open('./Python/Output/InversionAtOneImprecise.txt', 'w') as f:
            f.write('Taylor Uncertainty\tInput Uncertainty\tValue\tUncertainty\n')
            for dx in (0.05, 0.1, 0.15, 0.2):
                for i in range(0, 51, 2):
                    dy = i / 100
                    sTaylor = [VarDbl(1 if (i % 2) == 0 else -1, dy) for i in range(taylor.momentum.maxOrder)]
                    res = taylor.taylor1d(VarDbl(1, dx), 'Imprecise_Coeff', sTaylor, True, True)
                    f.write(f'{dy}\t{dx}\t{res.value()}\t{res.uncertainty()}\n')




if __name__ == '__main__':
    unittest.main()