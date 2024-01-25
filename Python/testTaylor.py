import math
import os
import numpy as np
import unittest
import sys

from taylor import Taylor, LossUncertaintyException
from varDbl import VarDbl, UncertaintyException, validate


taylor = Taylor()

HIST_RANGE = 3
HIST_DIVIDES = 5
HIST_BINS = 3 * 2 * HIST_DIVIDES
SAMPLES = 10000

def writeTest(self, test:str, func, npfun, sX:tuple[float],
              sDev=[0.2] + [math.pow(10,-err) for err in range(1, 17)],
              ignoreAssertError:bool=False ):
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
                    print(f'Ignore {test}({x}+/-{dev}) due to AssertionError {ex}')
                    continue
                except LossUncertaintyException as ex:
                    print(f'Ignore {test}({x}+/-{dev}) due to LossUncertaintyException {ex}')
                    continue
                except BaseException as ex:
                    print(f'Ignore {test}({x}+/-{dev}) due to BaseException {ex}')
                    continue
                
                for NoiseType, sNoise in ssNoise.items():
                    try:
                        sError = npfun(x, sNoise) - base
                    except Exception as ex:
                        print(f'Ignore {test}({x}+/-{dev}) due to numpy Exception {ex}')
                        continue
                    except BaseException:
                        print(f'Ignore {test}({x}+/-{dev}) due to numpy BaseException {ex}')
                        continue
                    if not math.isfinite(np.min(sError)):
                        print(f'Ignore {test}({x}+/-{dev}) due to numpy nan')
                        continue
                    valDev = np.std(sError)
                    sNorm = sError / var.uncertainty()
                    f.write(f"{NoiseType}\t{dev}\t{x}\t{base}\t{np.std(sNorm)}\t{np.min(sNorm)}\t{np.max(sNorm)}")
                    f.write(f"\t{valDev}\t{var.uncertainty()}\t{np.mean(sError)}\t{var.value() - base}")
                    hist, bin_edges = np.histogram(sNorm, bins=HIST_BINS, range=(-HIST_RANGE, +HIST_RANGE), density=True)
                    for i in range(HIST_BINS):
                        f.write(f'\t{hist[i]}')
                    f.write('\n')


class TestExp (unittest.TestCase):
    sTaylor = taylor.exp()

    def testCoeff(self):
        validate(self, TestExp.sTaylor[1], 1/1)
        validate(self, TestExp.sTaylor[2], 1/2)
        validate(self, TestExp.sTaylor[3], 1/6)
        validate(self, TestExp.sTaylor[4], 1/24)
        validate(self, TestExp.sTaylor[5], 1/120)

    def validate(self, exp, uncertainty, exception=None, 
                 valueDelta=5e-7, uncertaintyDelta=1e-6) -> VarDbl:
        try:
            s1dTaylor = TestExp.sTaylor
            s1dTaylor[0] = math.exp(exp)
            var = VarDbl(exp, uncertainty)
            res = taylor.taylor1d(var, "exp", s1dTaylor, False, True)
            prec = res * math.exp(-exp)
            if exception:
                self.fail(f'precision {prec} of exp({var}) = {res} not throw {exception}')
            self.assertAlmostEqual(prec.value(),
                    1 + math.pow(uncertainty, 2)/2 + math.pow(uncertainty, 4)/8 +\
                        math.pow(uncertainty, 6)/48 + math.pow(uncertainty, 8)/384,
                    delta=valueDelta)
            self.assertAlmostEqual(prec.variance(),
                    math.pow(uncertainty, 2) + math.pow(uncertainty, 4)*3/2 +\
                        math.pow(uncertainty, 6)*7/6 + math.pow(uncertainty, 8)*5/8,
                    delta=uncertaintyDelta)
            return res
        except BaseException as ex:
            if (exception is None) or (type(ex) != exception): 
                raise ex
            return res

    def testException(self):
        self.assertAlmostEqual(709.78, math.log(sys.float_info.max), delta=0.1)

        self.validate(-392, 0.1, UncertaintyException)
        self.validate(-390, 0.1, AssertionError)
        self.validate(-200, 0.1)

        self.validate(196, 0.1, UncertaintyException)
        self.validate(194, 0.1, AssertionError)
        self.validate(100, 0.1)
 
        self.validate(-392, 0.2, UncertaintyException)
        self.validate(-390, 0.2, AssertionError)
        self.validate(-200, 0.2, valueDelta=1e-3)

        self.validate(196, 0.2, UncertaintyException)
        self.validate(194, 0.2, AssertionError)
        self.validate(100, 0.2, valueDelta=1e-3)
 
        self.validate(-392, 1, UncertaintyException)
        self.validate(-390, 1, AssertionError)
        self.validate(-200, 1, valueDelta=1e-3, uncertaintyDelta=0.5)

        self.validate(196, 1, UncertaintyException)
        self.validate(194, 1, AssertionError)
        self.validate(100, 1, valueDelta=1e-3, uncertaintyDelta=0.5)

    @staticmethod
    def func(x):
        return math.exp(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.exp(sNoise + x)
 
    def test(self):
        writeTest(self, 'exp', TestExp.func, TestExp.npfunc,
                  (-100, -50, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 50, 100))





class TestLog (unittest.TestCase):
    sTaylor = taylor.log()

    def testCoeff(self):
        validate(self, TestLog.sTaylor[1], 1/1)
        validate(self, TestLog.sTaylor[2], -1/2)
        validate(self, TestLog.sTaylor[3], 1/3)
        validate(self, TestLog.sTaylor[4], -1/4)
        validate(self, TestLog.sTaylor[5], 1/5)

    def validate(self, x, uncertainty, exception=None, 
                 valueDelta=2e-3, uncertaintyDelta=5e-3) -> VarDbl:
        try:
            s1dTaylor = TestLog.sTaylor
            s1dTaylor[0] = math.log(x)
            var = VarDbl(x, uncertainty)
            res = taylor.taylor1d(var, "log", s1dTaylor, True, False)
            precOut = res - math.log(x)
            if (exception is not None) and (exception != AssertionError):
                self.fail(f'precision {precOut} of log({var}) = {res} not throw {exception}')
            precIn = abs(uncertainty/x)
            self.assertAlmostEqual(precOut.value(),
                    - math.pow(precIn, 2)/2 - math.pow(precIn, 4)/4 -\
                        math.pow(precIn, 6)/6 - math.pow(precIn, 8)/8,
                    delta=valueDelta)
            self.assertAlmostEqual(precOut.variance(),
                    math.pow(precIn, 2) + math.pow(precIn, 4)*9/8 +\
                        math.pow(precIn, 6)*119/24 + math.pow(precIn, 8)*991/32,
                    delta=uncertaintyDelta)
            return res
        except BaseException as ex:
            if (exception is None) or (type(ex) != exception): 
                raise ex
            return res

    def testException(self):
        self.assertAlmostEqual(709.78, math.log(sys.float_info.max), delta=0.1)

        self.validate(1/4, 0.2, LossUncertaintyException)
        self.validate(1/2, 0.2, AssertionError)
        self.validate(1, 0.2)

        self.validate(1/8, 0.1, LossUncertaintyException)
        self.validate(1/4, 0.1, AssertionError)
        self.validate(1/2, 0.1)

        self.validate(1/64, 0.01, LossUncertaintyException)
        self.validate(1/32, 0.01, AssertionError)
        self.validate(1/16, 0.01)

    @staticmethod
    def func(x):
        return math.log(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.log(sNoise + x)
 
    def test(self):
        writeTest(self, 'log', TestLog.func, TestLog.npfunc,
                  (1/32, 1/20, 1/16, 0.1, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                  ignoreAssertError=True)
        

        


if __name__ == '__main__':
    unittest.main()