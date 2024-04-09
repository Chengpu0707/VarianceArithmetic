import math
import logging
import numpy as np
import os
import unittest
import sys

from taylor import Taylor, NotReliableException, NotMonotonicException, NotStableException
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
                except (NotMonotonicException, NotReliableException, NotStableException) as ex:
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

        self.validate(1/2, 0.2, NotMonotonicException)
        self.validate(3/4, 0.2, AssertionError)
        self.validate(1, 0.2)

        self.validate(1/4, 0.1, NotMonotonicException)
        self.validate(1/2, 0.1, AssertionError)
        self.validate(1, 0.1)

        self.validate(1/64, 0.01, NotMonotonicException)
        self.validate(1/32, 0.01, NotMonotonicException)
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
                  sDev = [math.pow(10,-err) for err in range(1, 17)] + [0, 0.2, 0.5])



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
        try:
            self.validate(0.5, 0.23,    valueDelta=7e-3, varianceDelta=2e-2)    
        except NotMonotonicException:
            pass    
        try:
            self.validate(-1, 0.2,      valueDelta=5e-3, varianceDelta=7e-3) 
        except NotMonotonicException:
            pass
        self.validate(-1, 0.195,    valueDelta=5e-3, varianceDelta=7e-3)  
        try:
            self.validate(-1.5, 0.205,  valueDelta=2e-2, varianceDelta=0.6) 
        except NotMonotonicException:
            pass
        try:
            self.validate(-1.5, 0.2,    valueDelta=2e-2, varianceDelta=8e-2) 
        except NotMonotonicException:
            pass
        self.validate(-1.5, 0.195,  valueDelta=2e-2, varianceDelta=3e-2)  
        try:
            self.validate(-2, 0.2,      valueDelta=3e-2, varianceDelta=1.7)  
        except NotMonotonicException:
            pass
        self.validate(-2, 0.194,    valueDelta=3e-2, varianceDelta=3e-1)  

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
    HEADER= 'Name\tParameter\tInput Value\tInput Uncertainty'\
            '\tOrder\tOrder Value\tOrder Variance'\
            '\tExpansion Value\tExpansion Variance'\
            '\tNew Value Value\tNew Value Uncertainty'\
            '\tNew Variance Value\tNew Variance Uncertainty\n'

    def dumpExpansion(self, fw, name:str, param:float, s1dTaylor:list[float]):
        unc = s1dTaylor[0].uncertainty()
        value = VarDbl(s1dTaylor[0].value())
        variance = VarDbl()
        var = VarDbl(s1dTaylor[0].variance())
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
            fw.write(f'{name}\t{param}\t{s1dTaylor[0].value()}\t{unc}'
                     f'\t{n}\t{value.value()}\t{variance.value()}'
                     f'\t{varn.value()}\t{varn.variance()}'
                     f'\t{newValue.value()}\t{newValue.variance()}'
                     f'\t{newVariance.value()}\t{newVariance.variance()}\n')
            fw.flush()
            value += newValue
            variance += newVariance
            varn *= var
            if varn.value() == 0:
                break
        return VarDbl(value.value(), variance.value() + value.variance(), True)
    
    def testPower(self):
        with open('./Python/Output/PowerExpansion.txt', 'w') as fw:
            fw.write(TestExpansion.HEADER)
            for i in range(-20, 0, 2):
                exp = i/10
                s1dTaylor = taylor.power(exp)
                for j in range (190, 201):
                    s1dTaylor[0] = VarDbl(1, j/1000)
                    res = self.dumpExpansion(fw, 'power', exp, s1dTaylor)
                    try:
                        dir = s1dTaylor[0] ** (exp)
                        self.assertAlmostEqual(res.value(), dir.value())
                        self.assertAlmostEqual(res.variance(), dir.variance())
                    except NotMonotonicException:
                        pass




class TestPolyNearOne (unittest.TestCase):

    def test_uncertainty_bias(self):
        '''
        1/(1 + |x|) starts to diverge at 0.06 noise
        '''
        for i in range(1, 8):
            noise = i/100
            for x in (-0.6, -0.7):
                x = VarDbl(x, noise)
                res = 1/(1 - x)
                try:
                    expd = taylor.polynominal(x, [1] * (taylor._momentum._maxOrder - 1))
                    if i == 7:
                        self.assertAlmostEqual(res.value(), expd.value(), places=3)
                        self.assertAlmostEqual(res.variance(), expd.variance(), places=1)
                    else:
                        self.assertAlmostEqual(res.value(), expd.value(), places=5)
                        self.assertAlmostEqual(res.variance(), expd.variance(), places=3)
                except AssertionError as ex:
                    raise ex

        for i in range(1, 6):
            noise = i/100
            for x in (0.6, 0.7):
                x = VarDbl(x, noise)
                res = 1/(1 - x)
                try:
                    expd = taylor.polynominal(x, [1] * (taylor._momentum._maxOrder - 1))
                    print(f'{x}: Taylor={expd}, precise={1/(1 - x)}')
                    if i == 6:
                        self.assertAlmostEqual(res.value(), expd.value(), places=3)
                        self.assertAlmostEqual(res.variance(), expd.variance(), places=1)
                    else:
                        self.assertAlmostEqual(res.value(), expd.value(), places=5)
                        self.assertAlmostEqual(res.variance(), expd.variance(), places=3)
                except AssertionError as ex:
                    if i == 7:
                        continue
                    raise ex


    @unittest.skip('Too slow')
    def test_precise(self):
        LOG_PATH = './Python/Output/PolyNearOne.log'
        if os.path.isfile(LOG_PATH):
            os.remove(LOG_PATH)
        logging.basicConfig(filename=LOG_PATH, encoding='utf-8', level=logging.DEBUG,
                            format='%(asctime)s:%(levelname)s:%(message)s')
        
        with open(f'./Python/Output/PolyNearOne.txt', 'w') as f:
            f.write('X\tOrder\tReminder\tValue Error\tUncertainty\tResult ULP'
                    '\tPower\tULP Power\tRouding Error\tAccumulated Rounding Error\n')
            
            def calc(x):
                final = 1/(1 - x)
                valSum = 1
                errSum = 0
                pow = 1
                for order in range(1, taylor._momentum._maxOrder - 1):
                    pow *= x
                    valSum += pow
                    try:
                        res = taylor.polynominal(VarDbl(x), [1] * (order + 1)) - final
                    except NotReliableException as ex:
                        logger.warning(f'At x={x} order={order} encounter NotReliableException={ex}')
                        break
                    except UncertaintyException as ex:
                        logger.warning(f'At x={x} order={order} encounter UncertaintyException={ex}')
                        break
                    self.assertAlmostEqual(valSum - final, res.value())
                    ulp = math.ulp(valSum)
                    ulpPow = ulp*round(pow/ulp)
                    err = pow - ulpPow
                    errSum += err
                    f.write(f'{x}\t{order}\t{valSum - final}\t{res.value()}\t{res.uncertainty()}\t{ulp}'
                            f'\t{pow}\t{ulpPow}\t{err}\t{errSum}\n')
                    f.flush()

            for i in range(50, 75, 5):
                calc(i / 100)
                calc(-i / 100)


    def calc_imprecise(self, x, noise, f=None):
        try:
            val = 1/(1 - x)
        except NotReliableException as ex:
            logger.warning(f'At x={x} noise={noise} encounter NotReliableException={ex}')
            return
        except UncertaintyException as ex:
            logger.warning(f'At x={x} noise={noise} encounter UncertaintyException={ex}')
            return
        
        stbUnc = False
        byInc = False
        byLSV = False
        terminated = False
        byErr = False
        byRem = False

        prevVal = 0
        prevUnc = 0
        prevStbUnc = False
        prevByInc = False
        prevByLSV = False
        prevTerminated = False
        prevByErr = False
        prevByRem = False
        for order in range(1, taylor._momentum._maxOrder):
            try:
                res = taylor.polynominal(x, [1] * order)
            except NotReliableException as ex:
                logger.warning(f'At x={x} noise={noise} order={order} encounter NotReliableException={ex}')
                break
            except UncertaintyException as ex:
                logger.warning(f'At x={x} noise={noise} order={order} encounter UncertaintyException={ex}')
                break
            err = res - val

            unc = res.uncertainty()
            incVal = res.value() - prevVal 
            incUnc = unc - prevUnc
            prevVal = res.value()
            prevUnc = unc

            unc *= Taylor.TAU
            uncVal = val.uncertainty()*Taylor.TAU
            stbUnc = abs(incUnc) < unc
            byInc = (abs(incVal) < unc)
            byLSV = (abs(incVal) <= math.ulp(res.value()))
            terminated = stbUnc and (byInc or byLSV)
            byErr = (abs(err.value()) < uncVal)
            byRem = (abs(err.value()) < math.ulp(val.value()))
            if f:
                f.write(f'{x.value()}\t{noise}\t{order}\t{res.value()}\t{res.uncertainty()}'
                        f'\t{math.ulp(res.value())}'
                        f'\t{incVal}\t{incUnc}\t{err.value()}\t{err.uncertainty()}'
                        f'\t{stbUnc != prevStbUnc}'
                        f'\t{byInc != prevByInc}'
                        f'\t{byLSV != prevByLSV}'
                        f'\t{terminated != prevTerminated}'
                        f'\t{byErr != prevByErr}'
                        f'\t{byRem != prevByRem}'
                        '\n')
                f.flush()
            prevStbUnc = stbUnc
            prevByInc = byInc
            prevByLSV = byLSV
            prevTerminated = terminated
            prevByErr = byErr
            prevByRem = byRem

            if stbUnc and byInc and byLSV and byErr and byRem:
                self.assertAlmostEqual(unc, val.uncertainty()*Taylor.TAU)
                break  
        return stbUnc, byErr, byInc, byLSV, byRem 
        
    def test_imprecise_test(self):
        res = 1/VarDbl(0.4, 0.06)
        self.assertAlmostEqual(res.value(), 2.5605577330563207)
        self.assertAlmostEqual(res.uncertainty(), 0.41519016618465204)
        try:
            res = 1/VarDbl(0.3, 0.06)
            self.assertAlmostEqual(res.value(), 3.4874748968790157)
            self.assertAlmostEqual(res.uncertainty(), 0.8333325278975393)
        except NotMonotonicException:
            return
        try:
            res = 1/VarDbl(0.3, 0.07)
            self.assertAlmostEqual(res.value(), 527.564157634998)
            self.assertAlmostEqual(res.uncertainty(), 452.0300821331497)
        except NotMonotonicException:
            return
        res = 1/VarDbl(0.3)
        self.assertAlmostEqual(res.value(), 3.3333333333333335)
        self.assertAlmostEqual(res.uncertainty(), 4.388015458239491e-16)
        self.assertAlmostEqual(math.ulp(1/3), 4.388015458239491e-16)

        self.assertTupleEqual(self.calc_imprecise(VarDbl(0.7), 0), 
                              (True, False, True, True, False))
        self.assertTupleEqual(self.calc_imprecise(VarDbl(0.7, 0.01), 0.01), 
                              (True, True, True, True, False))
        self.assertTupleEqual(self.calc_imprecise(VarDbl(0.7, 0.01), 0.04), 
                              (True, True, True, True, False))
        
        self.assertTupleEqual(self.calc_imprecise(VarDbl(-0.7), 0), 
                              (True, False, True, True, False))
        self.assertTupleEqual(self.calc_imprecise(VarDbl(-0.7, 0.01), 0.01), 
                              (True, True, True, True, True))
        self.assertTupleEqual(self.calc_imprecise(VarDbl(-0.7, 0.01), 0.04), 
                              (True, True, True, True, True))
    
    @unittest.skip('10 minutes slow')
    def test_imprecise(self):
        LOG_PATH = './Python/Output/UncertaintyNearOne.log'
        if os.path.isfile(LOG_PATH):
            os.remove(LOG_PATH)
        logging.basicConfig(filename=LOG_PATH, encoding='utf-8', level=logging.DEBUG,
                            format='%(asctime)s:%(levelname)s:%(message)s')
        
        with open(f'./Python/Output/UncertaintyNearOne.txt', 'w') as f:
            f.write('X\tInput Uncertainty\tOrder\tExpansion Value\tExpansion Uncertainty\tExpansion ULP'
                    '\tValue Increment\tUncertainty Increment\tError Value\tError Uncertainty'
                    '\tStable Uncertainty\tBy Increment\tBy LSV\tTerminated\tBy Error\tBy Reminder\n')

            for x in (0.5, 0.6, 0.7):
                for sign in (-1, 1):
                    self.calc_imprecise(VarDbl(sign*x), 0, f)
                    for i in range(1, 10):
                        noise = i/100
                        self.calc_imprecise(VarDbl(sign*x, noise), noise, f)
                        noise = i/1000
                        self.calc_imprecise(VarDbl(sign*x, noise), noise, f)



if __name__ == '__main__':
    unittest.main()