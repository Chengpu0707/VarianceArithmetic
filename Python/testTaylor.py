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

def dump_test(self, test:str, func, npfun, sX:tuple[float],
              sDev=[0.2] + [math.pow(10,-err) for err in range(1, 17)],
              ignoreAssertError:bool=True ):
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
                if var.uncertainty() <= 0:
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

    def test_taylor(self):
        validate(self, TestExp.sTaylor[1], 1/1, 1.2819751242557095e-16)
        validate(self, TestExp.sTaylor[2], 1/2, 9.06493303673679e-17)
        validate(self, TestExp.sTaylor[3], 1/6, 3.42026916246155e-17)
        validate(self, TestExp.sTaylor[4], 1/24, 1.008197910236068e-17)
        validate(self, TestExp.sTaylor[5], 1/120, 2.1240690246726458e-18)

    def validate(self, exp, uncertainty, exception=None, 
                 valueDelta=5e-7, varianceDelta=1e-6) -> VarDbl:
        try:
            s1dTaylor = TestExp.sTaylor
            s1dTaylor[0] = math.exp(exp)
            var = VarDbl(exp, uncertainty)
            res = taylor.taylor1d(var, "exp", s1dTaylor, False, True)
            prec = res * math.exp(-exp)
            if (exception is not None) and (exception != AssertionError):
                self.fail(f'precision {prec} of exp({var}) = {res} not throw {exception}')
            self.assertAlmostEqual(prec.value(),
                    1 + math.pow(uncertainty, 2)/2 + math.pow(uncertainty, 4)/8 +\
                        math.pow(uncertainty, 6)/48 + math.pow(uncertainty, 8)/384,
                    delta=valueDelta)
            self.assertAlmostEqual(prec.variance(),
                    math.pow(uncertainty, 2) + math.pow(uncertainty, 4)*3/2 +\
                        math.pow(uncertainty, 6)*7/6 + math.pow(uncertainty, 8)*5/8,
                    delta=varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res
        except BaseException as ex:
            if (exception is None) or (not isinstance(ex, exception)): 
                raise ex


    def test_exception(self):
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
        self.validate(-200, 1, valueDelta=1e-3, varianceDelta=0.5)

        self.validate(196, 1, UncertaintyException)
        self.validate(194, 1, AssertionError)
        self.validate(100, 1, valueDelta=1e-3, varianceDelta=0.5)

    @staticmethod
    def func(x):
        return math.exp(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.exp(sNoise + x)
 
    def test_dump(self):
        dump_test(self, 'exp', TestExp.func, TestExp.npfunc,
                  (-100, -50, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 50, 100),
                  sDev = [math.pow(10,-err) for err in range(0, 17)])





class TestLog (unittest.TestCase):
    sTaylor = taylor.log()

    def test_taylor(self):
        validate(self, TestLog.sTaylor[1], 1/1, 1.2819751242557095e-16)
        validate(self, TestLog.sTaylor[2], -1/2, 6.409875621278547e-17)
        validate(self, TestLog.sTaylor[3], 1/3, 3.2049378106392736e-17)
        validate(self, TestLog.sTaylor[4], -1/4, 3.2049378106392736e-17)
        validate(self, TestLog.sTaylor[5], 1/5, 1.6024689053196368e-17)

    def validate(self, x, uncertainty, exception=None, 
                 valueDelta=2e-3, varianceDelta=5e-3) -> VarDbl:
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
                    delta=varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res
        except BaseException as ex:
            if (exception is None) or (not isinstance(ex, exception)): 
                raise ex

    def test_exception(self):
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
 
    def test_dump(self):
        dump_test(self, 'log', TestLog.func, TestLog.npfunc,
                  (1/32, 1/20, 1/16, 0.1, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32),
                  ignoreAssertError=True)
        

class TestSin (unittest.TestCase):

    def test_taylor_0(self):
        sTaylor = taylor.sin(0)
        validate(self, sTaylor[0], 0)
        validate(self, sTaylor[1], 1/1, 0.8164965809277261*math.ulp(1))
        validate(self, sTaylor[2], 0)
        validate(self, sTaylor[3], -1/6, 1.4529663145135579*math.ulp(1/6))
        validate(self, sTaylor[4], 0)
        validate(self, sTaylor[5], 1/120, 1.3705905728986023*math.ulp(1/120))

    def test_taylor_half_pi(self):
        sTaylor = taylor.sin(math.pi/2)
        validate(self, sTaylor[0], 1, 0.57735026918962581*math.ulp(1))
        validate(self, sTaylor[1], 0, 1.059541960031773e-32, deltaValue=6.123233995736766e-17)
        validate(self, sTaylor[2], -1/2)
        validate(self, sTaylor[3], 0, 2.4068420546318608e-33, deltaValue=1.020538999289461e-17)
        validate(self, sTaylor[4], 1/24, 1.6442942874387492*math.ulp(1/24))
        validate(self, sTaylor[5], 0, 1.4294378989955562e-34, deltaValue=5.102694996447305e-19)
       
    def test_taylor_half_pi_negative(self):
        sTaylor = taylor.sin(-math.pi/2)
        validate(self, sTaylor[0], -1, 0.57735026918962581*math.ulp(1))
        validate(self, sTaylor[1], 0, 1.059541960031773e-32, deltaValue=6.123233995736766e-17)
        validate(self, sTaylor[2], 1/2)
        validate(self, sTaylor[3], 0, 2.4068420546318608e-33, deltaValue=1.020538999289461e-17)
        validate(self, sTaylor[4], -1/24, 1.6442942874387492*math.ulp(1/24))
        validate(self, sTaylor[5], 0, 1.4294378989955562e-34, deltaValue=5.102694996447305e-19)
       
    def validate(self, x, uncertainty, exception=None, 
                 valueDelta=5e-7, varianceDelta=1e-6) -> VarDbl:
        try:
            s1dTaylor = taylor.sin(x)
            s1dTaylor[0] = math.sin(x)
            var = VarDbl(x, uncertainty)
            res = taylor.taylor1d(var, "sin", s1dTaylor, False, False)
            prec = res - s1dTaylor[0]
            if (exception is not None) and (exception != AssertionError):
                self.fail(f'precision {prec} of sin({var}) = {res} not throw {exception}')
            self.assertAlmostEqual(prec.value(),
                    math.sin(x) * (-math.pow(uncertainty, 2)/2 + math.pow(uncertainty, 4)/8 +\
                        -math.pow(uncertainty, 6)/48 + math.pow(uncertainty, 8)/384),
                    delta=valueDelta)
            cos2 = math.cos(x) * math.cos(x)
            self.assertAlmostEqual(prec.variance(),
                    math.pow(uncertainty, 2) * cos2 - math.pow(uncertainty, 4)*(3/2 * cos2 - 1/2) +\
                        math.pow(uncertainty, 6)*(7/6 * cos2 - 1/2),
                    delta=varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res
        except BaseException as ex:
            if (exception is None) or (not isinstance(ex, exception)): 
                raise ex
        
    def test_exception(self):
        self.validate(0, 0)
        self.validate(0, 1e-3)
        self.validate(0, 1e-2)
        self.validate(0, 1e-1)

        self.validate(math.pi/2, 0)
        self.validate(math.pi/2, 1e-3)
        self.validate(math.pi/2, 1e-2)
        self.validate(math.pi/2, 1e-1)

        self.validate(-math.pi/2, 0)
        self.validate(-math.pi/2, 1e-3)
        self.validate(-math.pi/2, 1e-2)
        self.validate(-math.pi/2, 1e-1)

    @staticmethod
    def func(x):
        return math.sin(x)

    @staticmethod
    def npfunc(x, sNoise):
        return np.sin(sNoise + x)
 
    def test_dump(self):
        dump_test(self, 'sin', TestSin.func, TestSin.npfunc,
                  np.array([i/16 for i in range(-16, 17)]) * math.pi,
                  sDev = [math.pow(10,-err) for err in range(0, 17)])



class TestPower (unittest.TestCase):

    def test_taylor_0(self):
        sTaylor = taylor.power(0)
        for coeff in sTaylor:
            validate(self, coeff, 0)

    def test_taylor_1(self):
        sTaylor = taylor.power(1)
        validate(self, sTaylor[1], 1, 0)
        for coeff in sTaylor[2:]:
            validate(self, coeff, 0, 0)

    def test_taylor_2(self):
        sTaylor = taylor.power(2)
        validate(self, sTaylor[1], 2, 0)
        validate(self, sTaylor[2], 1, 1.2819751242557095e-16)
        for coeff in sTaylor[3:]:
            validate(self, coeff, 0, 0)

    def test_taylor_sqrt(self):
        sTaylor = taylor.power(0.5)
        validate(self, sTaylor[1],  1/2,   6.409875621278547e-17)
        validate(self, sTaylor[2], -1/8,   2.2662332591841976e-17)
        validate(self, sTaylor[3],  1/16,  1.3877787807814457e-17)
        validate(self, sTaylor[4], -5/128, 9.554111923975634e-18)
        validate(self, sTaylor[5],  7/256, 7.141219782764964e-18)

    def test_taylor_inv_sqrt(self):
        sTaylor = taylor.power(-0.5)
        validate(self, sTaylor[1],  -1/2,   6.409875621278547e-17)
        validate(self, sTaylor[2],   3/8,   5.777783805466599e-17)
        validate(self, sTaylor[3],  -5/16,  5.381475625187719e-17)
        validate(self, sTaylor[4],  35/128, 5.117134779608756e-17)
        validate(self, sTaylor[5], -63/256, 4.9276631313251286e-17)

    def test_taylor_inv_1(self):
        sTaylor = taylor.power(-1)
        validate(self, sTaylor[1], -1, 0)
        validate(self, sTaylor[2],  1, 1.2819751242557095e-16)
        validate(self, sTaylor[3], -1, 1.812986607347358e-16)
        validate(self, sTaylor[4],  1, 2.220446049250313e-16)
        validate(self, sTaylor[5], -1, 2.563950248511419e-16)

    def test_taylor_inv_2(self):
        sTaylor = taylor.power(-2)
        validate(self, sTaylor[1], -2, 0)
        validate(self, sTaylor[2],  3, 2.563950248511419e-16)
        validate(self, sTaylor[3], -4, 5.145674902128044e-16)
        validate(self, sTaylor[4],  5, 8.226007047307469e-16)
        validate(self, sTaylor[5], -6, 1.1769760485126625e-15)

    def validate(self, exp, uncertainty, exception=None, 
                 valueDelta=1e-3, varianceDelta=5e-3) -> VarDbl:
        try:
            s1dTaylor = taylor.power(exp)
            s1dTaylor[0] = VarDbl(1, 0)
            var = VarDbl(1, uncertainty)
            res = taylor.taylor1d(var, "pow", s1dTaylor, True, True)
            prec = res
            if (exception is not None) and (exception != AssertionError):
                self.fail(f'precision {prec} of exp({var}) = {res} not throw {exception}')
            self.assertAlmostEqual(prec.value(),
                    1 + math.pow(uncertainty, 2) * exp*(exp-1)/2 + \
                        math.pow(uncertainty, 4) * exp*(exp-1)*(exp-2)*(exp-3)/24 +\
                        math.pow(uncertainty, 6) * exp*(exp-1)*(exp-2)*(exp-3)*(exp-4)*(exp-5)/720,
                    delta=valueDelta)
            self.assertAlmostEqual(prec.variance(),
                    math.pow(uncertainty, 2) * exp*exp + \
                    math.pow(uncertainty, 4) * exp*exp*(exp-1)*(exp-5/3)*3/2 +\
                    math.pow(uncertainty, 6) * exp*exp*(exp-1)*(exp-2)*(exp-2)*(exp-16/7)*7/6,
                    delta=varianceDelta)
            return res
        except AssertionError as ex:
            if (exception is None) or (not isinstance(ex, AssertionError)): 
                raise ex
            else:
                return res
        except BaseException as ex:
            if (exception is None) or (not isinstance(ex, exception)): 
                raise ex
            
    def test_near_2(self):
        two = self.validate(2, 0.2)
        validate(self, two, 1.039999427727274, 0.4039766494889215)
        lower = self.validate(2 - 1e-9, 0.2)
        validate(self, lower, two.value(), two.uncertainty(), deltaValue=6e-11, deltaUncertainty=3e-10)
        upper = self.validate(2 + 1e-9, 0.2)
        validate(self, upper, two.value(), two.uncertainty(), deltaValue=6e-11, deltaUncertainty=3e-10)


    def test_exception(self):
        self.validate(0, 0) 
        self.validate(0, 0.2) 
        self.validate(1e-6, 0)    
        self.validate(1e-3, 0)   
        self.validate(0.1, 0)   

        self.validate(1, 0.2)  
        self.validate(2, 0.2)    
        self.validate(0.5, 0.2)    
        self.validate(0.5, 0.205)    
        self.validate(0.5, 0.23,    valueDelta=7e-3, varianceDelta=2e-2)    
        self.validate(-1, 0.2,      valueDelta=5e-3, varianceDelta=7e-3) 
        self.validate(-1, 0.195,    valueDelta=5e-3, varianceDelta=7e-3)  
        self.validate(-1.5, 0.205,  valueDelta=2e-2, varianceDelta=0.6) 
        self.validate(-1.5, 0.2,    valueDelta=2e-2, varianceDelta=8e-2) 
        self.validate(-1.5, 0.195,  valueDelta=2e-2, varianceDelta=3e-2)  
        self.validate(-2, 0.2,      valueDelta=3e-2, varianceDelta=1.7)  
        self.validate(-2, 0.195,    valueDelta=3e-2, varianceDelta=3e-1)  

    @staticmethod
    def func(x):
        return 1.0

    @staticmethod
    def npfunc(x, sNoise):
        return np.power(sNoise + 1, x)
 
    def test_dump(self):
        dump_test(self, 'pow', TestPower.func, TestPower.npfunc,
                  np.array([i/10 for i in range(-20, 31, 2)] +[-0.1, 0.1, -0.01, 0.01, -1e-3, 1e-6]),
                  sDev = [0.2, 0.195] + [math.pow(10,-err) for err in range(1, 17)])


if __name__ == '__main__':
    unittest.main()