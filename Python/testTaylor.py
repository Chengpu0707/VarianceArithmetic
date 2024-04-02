import math
import logging
import numpy as np
import os
import unittest
import sys

from taylor import Taylor, LossUncertaintyException
from varDbl import VarDbl, UncertaintyException, validate

logger = logging.getLogger(__name__)

taylor = Taylor()

HIST_RANGE = 3
HIST_DIVIDES = 5
HIST_BINS = 3 * 2 * HIST_DIVIDES
SAMPLES = 10000

def dump_test(self, test:str, func, npfun, sX:tuple[float],
              sDev=[0.2] + [math.pow(10,-err) for err in range(1, 17)] + [0],
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
            var = VarDbl(exp, uncertainty)
            res = VarDbl.exp(var)
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
            var = VarDbl(x, uncertainty)
            res = VarDbl.log(var)
            precOut = res - math.log(var.value())
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

        self.validate(1/2, 0.2, LossUncertaintyException)
        self.validate(3/4, 0.2, AssertionError)
        self.validate(1, 0.2)

        self.validate(1/4, 0.1, LossUncertaintyException)
        self.validate(1/2, 0.1, AssertionError)
        self.validate(1, 0.1)

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
            var = VarDbl(x, uncertainty)
            res = VarDbl.sin(var)
            prec = res - math.sin(var.value())
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
                  sDev = [math.pow(10,-err) for err in range(0, 17)] + [0])



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
        validate(self, sTaylor[2], 1, VarDbl.ulp(1))
        for coeff in sTaylor[3:]:
            validate(self, coeff, 0, 0)

    def test_taylor_3(self):
        sTaylor = taylor.power(3)
        validate(self, sTaylor[1], 3, 0)
        validate(self, sTaylor[2], 3, VarDbl.ulp(3)*1.5)
        validate(self, sTaylor[3], 1, VarDbl.ulp(1)*1.25, deltaValue=2.23e-16)
        for coeff in sTaylor[4:]:
            validate(self, coeff, 0, 0)

    def test_taylor_4(self):
        sTaylor = taylor.power(4)
        validate(self, sTaylor[1], 4, 0)
        validate(self, sTaylor[2], 6, VarDbl.ulp(6))
        validate(self, sTaylor[3], 4, VarDbl.ulp(4)*1.003466214899358)
        validate(self, sTaylor[4], 1, VarDbl.ulp(1)*1.4166666666666667)
        for coeff in sTaylor[5:]:
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
            var = VarDbl(1, uncertainty)
            res = var ** exp
            prec = res - VarDbl(1, 0)
            if (exception is not None) and (exception != AssertionError):
                self.fail(f'precision {prec} of exp({var}) = {res} not throw {exception}')
            self.assertAlmostEqual(prec.value(),
                    math.pow(uncertainty, 2) * exp*(exp-1)/2 + \
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


class TestPolynominial (unittest.TestCase):

    def test_poly_0(self):
        res = taylor.polynominal(VarDbl(), (1,))
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

        res = taylor.polynominal(VarDbl(1), (1,))
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

        res = taylor.polynominal(VarDbl(1), (-2,))
        self.assertEqual(res.value(), -2)
        self.assertEqual(res.uncertainty(), 0)

        res = taylor.polynominal(VarDbl(2), (VarDbl(1.0),))
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), VarDbl.ulp(1))

    def test_poly_1(self):
        res = taylor.polynominal(VarDbl(0, 1/8), (0,1))
        self.assertEqual(res.value(), 0)
        self.assertAlmostEqual(res.uncertainty(), 1/8, delta=1e-6)

        res = taylor.polynominal(VarDbl(0, 1/8), (1,1))
        self.assertEqual(res.value(), 1)
        self.assertAlmostEqual(res.uncertainty(), 1/8, delta=1e-6)

        res = taylor.polynominal(VarDbl(-2, 1/8), (1,1))
        self.assertEqual(res.value(), -1)
        self.assertAlmostEqual(res.uncertainty(), 1/8, delta=1e-6)

    def test_poly_2(self):
        for value in (0, 1, -1, 2, -2, 0.25, -0.25):
            try:
                res = taylor.polynominal(VarDbl(value, 0.5), (0,0,1))
                res2 = VarDbl(value, 0.5)**2
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            try:
                res = taylor.polynominal(VarDbl(value, 0.5), (1,2,1))
                res2 = VarDbl(1 + value, 0.5)**2
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            try:
                res = taylor.polynominal(VarDbl(value, 0.5), (1,-2,1))
                res2 = VarDbl(1 - value, 0.5)**2
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex

    def test_poly_3(self):
        for value in (0, 1, -1, 2, -2, 0.25, -0.25):
            try:
                res = taylor.polynominal(VarDbl(value, 0.5), (0,0,0,1))
                res2 = VarDbl(value, 0.5)**3
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            try:
                res = taylor.polynominal(VarDbl(value, 0.5), (1,3,3,1))
                res2 = VarDbl(1 + value, 0.5)**3
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            try:
                res = taylor.polynominal(VarDbl(value, 0.5), (1,-3,3,-1))
                res2 = VarDbl(1 - value, 0.5)**3
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            
    def varify_near_one(self, x, order):
        sPlus = [1, -1] * ((order // 2) + 1)
        if len(sPlus) > (order + 1):
            sPlus.pop(-1)
        self.assertEqual(len(sPlus), order + 1)
        for i in range(order):
            self.assertEqual(sPlus[i], -1 if i % 2 else 1)
        sMinus = [1] * (order + 1)
        
        plus = 0
        minus = 0
        pow = 1
        for i in range(order + 1):
            plus += sPlus[i] * pow
            minus += sMinus[i] * pow
            pow *= x
            if pow < max(math.ulp(plus), math.ulp(minus)):
                break
        resPlus = taylor.polynominal(VarDbl(x), sPlus)
        self.assertAlmostEqual(resPlus.value(), plus)
        resMinus = taylor.polynominal(VarDbl(x), sMinus)
        self.assertAlmostEqual(resMinus.value(), minus)
        return resPlus, resMinus

            
    def test_near_one(self):
        self.varify_near_one(0.7, 100)

        with open(f'./Python/Output/PolyNearOne.txt', 'w') as f:
            f.write('X\tOrder\tValue Error\tUncertainty'
                    '\tPower\tReminder\tRouding Error\tAccumulated Rounding Error'
                    '\tResult ULP\tULP Power\tDecrease\tULP Reminder\n')
            
            def calc(x):
                plusVal = 1/(1 + x)
                minusVal = 1/(1 - x)

                plus = 1
                minus = 1
                pow = x
                remPlus = 1 - plusVal
                remMinus = 1 - minusVal
                acmPlus = 0
                acmMinus = 0
                for order in range(1, taylor._momentum._maxOrder - 1):
                    plus += (-1 if (order % 2) else 1) * pow
                    minus += pow
                    plusRes, minusRes = self.varify_near_one(x, order)

                    res = plusRes - plusVal
                    rem = plus - plusVal
                    rounding = rem - remPlus
                    acmPlus += (pow - abs(rounding)) if (order % 2) else (abs(rounding) - pow)
                    ulp = math.ulp(plus)
                    f.write(f'{+x}\t{order}\t{res.value()}\t{res.uncertainty()}'
                            f'\t{pow}\t{rem}\t{rounding}\t{acmPlus}'
                            f'\t{ulp}\t{ulp*round(pow/ulp)}\t{rem - remPlus}\t{ulp*round(rem/ulp)}'
                            '\n')
                    remPlus = rem
                    
                    res = minusRes - minusVal
                    rem = minus - minusVal
                    rounding = remMinus - rem
                    acmMinus += abs(rounding) - pow
                    ulp = math.ulp(minus)
                    f.write(f'{-x}\t{order}\t{res.value()}\t{res.uncertainty()}'
                            f'\t{pow}\t{rem}\t{rounding}\t{acmMinus}'
                            f'\t{ulp}\t{ulp*round(pow/ulp)}\t{rem - remPlus}\t{ulp*round(rem/ulp)}'
                            '\n')
                    remMinus = rem                    

                    f.flush()
                    pow *= x

            for i in range(50, 75, 5):
                calc(i / 100)

    #@unittest.skip('Too slow')
    def test_near_one_uncertainty(self):
        LOG_PATH = './Python/Output/UncertaintyNearOne.log'
        if os.path.isfile(LOG_PATH):
            os.remove(LOG_PATH)
        logging.basicConfig(filename=LOG_PATH, encoding='utf-8', level=logging.DEBUG,
                            format='%(asctime)s:%(levelname)s:%(message)s')
        
        with open(f'./Python/Output/UncertaintyNearOne.txt', 'w') as f:
            f.write('X\tInput Uncertainty\tOrder\tExpansion Value\tExpansion Uncertainty'
                    '\tExpansion ULP\tError Value\tError Uncertainty\n')

            def calc(x, noise):
                try:
                    val = 1/(1 - x)
                except LossUncertaintyException as ex:
                    logger.warning(f'At x={x} noise={noise} encounter LossUncertaintyException={ex}')
                    return
                except UncertaintyException as ex:
                    logger.warning(f'At x={x} noise={noise} encounter UncertaintyException={ex}')
                    return
                
                for order in range(1, taylor._momentum._maxOrder - 1):
                    try:
                        res = taylor.polynominal(x, [1] * (order + 1))
                    except LossUncertaintyException as ex:
                        logger.warning(f'At x={x} noise={noise} order={order} encounter LossUncertaintyException={ex}')
                        break
                    except UncertaintyException as ex:
                        logger.warning(f'At x={x} noise={noise} order={order} encounter UncertaintyException={ex}')
                        break
                    err = res - val
                    f.write(f'{x.value()}\t{noise}\t{order}\t{res.value()}\t{res.uncertainty()}'
                            f'\t{math.ulp(res.value())}\t{err.value()}\t{err.uncertainty()}\n')
                    f.flush()
                    if abs(err.value()) <= math.ulp(res.value()):
                          break                           

            for x in (0.5, 0.6, 0.7):
                for sign in (-1, 1):
                    calc(VarDbl(sign*x), 0)
                    for noise in (1e-3, 1e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2):
                        calc(VarDbl(sign*x, noise), noise)

    def test_near_one_uncertainty_bias(self):
        '''
        1/(1 - 6.000000e-01~1.000e-02)==2.501565e+00~6.266e-02 vs 2.5
        1/(1 - 6.000000e-01~5.000e-02)==2.541054e+00~3.344e-01 vs 2.5
        1/(1 - 6.000000e-01~7.000e-02)==2.585001e+00~5.070e-01 vs 2.5
        1/(1 - 6.000000e-01~9.000e-02)==8.228626e+00~3.997e+01 vs 2.5
        1/(1 - 6.000000e-01~1.000e-01)==1.478180e+06~2.099e+04 vs 2.5
        '''
        for noise in (0.01, 0.05, 0.09, 0.1):
            for x in (0.6,):
                res = 1/(1 - VarDbl(x, noise))
                try:
                    self.assertAlmostEqual(res.value(), 1/(1 - x))
                    self.fail(f'x={x}: {res}=={1/(1 - x)}')
                except AssertionError:
                    print(f'1/(1 - {VarDbl(x, noise)})=={res} vs {1/(1 - x)}')


class TestLibError (unittest.TestCase):
    
    def testExpLog(self):
        with open(f'./Python/Output/ExpLogError.txt', 'w') as f:
            f.write('X\tVarDbl Error\tVarDbl Uncertainty\tLib Error\tLib Uncertainty\n')
            for i in range(100):
                x = i / 50
                var = VarDbl.log(VarDbl.exp(VarDbl(float(-x)))) + x
                lib = VarDbl(math.log(math.exp(-x))) + x
                f.write(f'{-x}\t{var.value()}\t{var.uncertainty()}\t{lib.value()}\t{lib.uncertainty()}\n')
                var = VarDbl.log(VarDbl.exp(VarDbl(float(x)))) - x
                lib = VarDbl(math.log(math.exp(x))) - x
                f.write(f'{x}\t{var.value()}\t{var.uncertainty()}\t{lib.value()}\t{lib.uncertainty()}\n')

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
                        var = (x **(exp))**(1/exp) - x
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


class TestExpansion (unittest.TestCase):
    HEADER= 'Name\tX\tOrder\tValue\tVariance'\
            '\tExpansion Value\tExpansion Variance'\
            '\tNew Value Value\tNew Value Uncertainty'\
            '\tNew Variance Value\tNew Variance Uncertainty\n'

    def dumpExpansion(self, fw, name:str, x: float, value:VarDbl, s1dTaylor:list[float]):
        variance = VarDbl()
        var = VarDbl(VarDbl(x).variance())
        varn = VarDbl(var)
        for n in range(2, taylor._momentum._maxOrder, 2):
            newValue = s1dTaylor[n] * taylor._momentum.factor(n) * varn
            newVariance = VarDbl()
            for j in range(1, n):
                newVariance += s1dTaylor[j] * s1dTaylor[n - j] * taylor._momentum.factor(n) * varn
            for j in range(2, n, 2):
                newVariance -= s1dTaylor[j] * taylor._momentum.factor(j) * \
                            s1dTaylor[n - j] * taylor._momentum.factor(n - j) * \
                            varn
            fw.write(f'{name}\t{x}\t{n}\t{value.value()}\t{variance.value()}')
            fw.write(f'\t{varn.value()}\t{varn.variance()}')
            fw.write(f'\t{newValue.value()}\t{newValue.variance()}')
            fw.write(f'\t{newVariance.value()}\t{newVariance.variance()}\n')
            value += newValue
            variance += newVariance
            varn *= var
            if varn.value() == 0:
                break
        return VarDbl(value.value(), variance.value() + value.variance(), True)


    def testSin(self):
        '''
        A reproduction of the logic in Tayor.taylor1d()
        '''
        with open('./Python/Output/Sin_Taylor.txt', 'w') as fw:
            fw.write(TestExpansion.HEADER)
            for i in range(1, 10):
                x = math.pi/4 *i
                value = VarDbl(math.sin(x))
                s1dTaylor = taylor.sin(value.value())
                res = self.dumpExpansion(fw, 'sin', x, value, s1dTaylor)
                self.assertAlmostEqual(math.sin(x), res.value(), delta=res.uncertainty())

    def testInversionNearZero(self):
        '''
        A reproduction of the logic in Tayor.taylor1d()
        '''
        with open('./Python/Output/NearOne_Taylor.txt', 'w') as fw:
            fw.write(TestExpansion.HEADER)
            for i in range(1, 16):
                x = math.pow(10, -i)
                s1dTaylor = taylor.power(-1)
                res = self.dumpExpansion(fw, '1/x', x, VarDbl(1, 0), s1dTaylor)
                self.assertAlmostEqual(1, res.value(), delta=res.uncertainty())


if __name__ == '__main__':
    unittest.main()