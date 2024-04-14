from contextlib import AbstractContextManager
import math
import pickle
import unittest
import random
import sys

from histo import Stat, Histo
from varDbl import VarDbl, ValueException, UncertaintyException, validate
from taylor import NotReliableException, NotMonotonicException


class TestInit (unittest.TestCase):

    def testInt(self):
        validate(self, VarDbl(0), 0, 0)
        validate(self, VarDbl(1), 1, 0)
        validate(self, VarDbl(-1), -1, 0)


    def testLargeInt(self):
        validate(self, VarDbl((1 << 53)), (1 << 53), 0)
        validate(self, VarDbl(-(1 << 53)), -(1 << 53), 0)

        # lost resolution
        f = float((1 << 53) + 1)
        validate(self, VarDbl((1 << 53) + 1), f, VarDbl.ulp(f))
        # no lost of resolution
        f = float(1 << 54)
        validate(self, VarDbl(1 << 54), f, 0)

    def testInit(self):
        validate(self, VarDbl(-1, 0), -1, 0)
        validate(self, VarDbl(-1, 1), -1, 1)
        validate(self, VarDbl(-1, -1), -1, 1)

    def failValueException(self, value):
        try:
            VarDbl(value)
            self.fail(f'Init ValDbl with {value}')
        except ValueException as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)

    def testValueException(self):
        self.failValueException(float('nan'))
        self.failValueException(float('inf'))

    def failUncertaintyException(self, value, uncertainty):
        try:
            VarDbl(value, uncertainty)
            self.fail(f'Init ValDbl with {value}~{uncertainty}')
        except UncertaintyException as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)

    def testUncertaintyException(self):
        self.failUncertaintyException(0, float('nan'))
        self.failUncertaintyException(0, float('inf'))
            
    def testUncertaintyRange(self):
        maxU = math.sqrt(sys.float_info.max)
        self.failUncertaintyException(0, maxU + VarDbl.ulp(maxU))
        validate(self, VarDbl(0, maxU), 0, maxU) 
        
        minU = math.sqrt(VarDbl.ulp(sys.float_info.min))
        validate(self, VarDbl(0, minU), 0, minU) 
        validate(self, VarDbl(0, minU*0.5), 0, 0) 


class TestRepresentation (unittest.TestCase):
    def testStr(self):
        v = VarDbl(-math.sqrt(2), math.sqrt(2))
        self.assertEqual('-1.414214e+00~1.414e+00', str(v))

    def testRepr(self):
        v = VarDbl(-math.sqrt(2), math.sqrt(2))
        with open('./Java/Output/data.pickle', 'wb') as f:
            pickle.dump(v, f, pickle.HIGHEST_PROTOCOL)
        with open('./Java/Output/data.pickle', 'rb') as f:
            vr = pickle.load(f)
        validate(self, vr, v.value(), v.uncertainty())

    def testBool(self):
        self.assertTrue(VarDbl(1))
        self.assertTrue(VarDbl(0, 1))
        self.assertFalse(VarDbl())
 

class TestAddSub (unittest.TestCase):
    def testNeg(self):
        validate(self, -VarDbl(1, 2), -1, 2) 

    def testAddVarDbl(self):
        v1 = VarDbl(math.sqrt(2), math.sqrt(2))
        v2 = VarDbl(-math.sqrt(2), math.sqrt(2))
        validate(self, v1 + v2, 0, 2) 
        v1 += v2
        validate(self, v1, 0, 2) 

    def testSubVarDbl(self):
        v1 = VarDbl(math.sqrt(2), math.sqrt(2))
        v2 = VarDbl(-math.sqrt(2), math.sqrt(2))
        validate(self, v2 - v1, -2*math.sqrt(2), 2) 
        v2 -= v1
        validate(self, v2, -2*math.sqrt(2), 2) 

    def testAddInt(self):
        v1 = VarDbl(1, math.sqrt(2))
        validate(self, v1 + 2, 3, math.sqrt(2)) 
        validate(self, 2 + v1, 3, math.sqrt(2)) 
        v1 += 2
        validate(self, v1, 3, math.sqrt(2)) 

    def testSubInt(self):
        v1 = VarDbl(1, math.sqrt(2))
        validate(self, v1 - 2, -1, math.sqrt(2)) 
        validate(self, 2 - v1, 1, math.sqrt(2)) 
        v1 -= 2
        validate(self, v1, -1, math.sqrt(2)) 

    def testAddFloat(self):
        v1 = VarDbl(1, math.sqrt(2))
        validate(self, v1 + 2.0, 3, math.sqrt(2)) 
        validate(self, 2.0 + v1, 3, math.sqrt(2)) 
        v1 += 2.0
        validate(self, v1, 3, math.sqrt(2)) 

        v2 = VarDbl(1.0)
        uncertainty = VarDbl.ulp(1) * math.sqrt(5)
        validate(self, v2 + 2.0, 3, uncertainty) 
        validate(self, 2.0 + v2, 3, uncertainty) 
        v2 += 2.0
        validate(self, v2, 3, uncertainty) 

    def testSubFloat(self):
        v1 = VarDbl(1, math.sqrt(2))
        validate(self, v1 - 2.0, -1, math.sqrt(2)) 
        v = 2.0 - v1
        validate(self, 2.0 - v1, 1, math.sqrt(2)) 
        v1 -= 2.0
        validate(self, v1, -1, math.sqrt(2)) 

        v2 = VarDbl(1.0)
        uncertainty = VarDbl.ulp(1) * math.sqrt(5)
        validate(self, v2 - 2.0, -1,  uncertainty) 
        validate(self, 2.0 - v2, 1,  uncertainty) 
        v2 -= 2.0
        validate(self, v2, -1,  uncertainty) 

    def testAddSubException(self):
        maxV = sys.float_info.max
        maxU = math.sqrt(sys.float_info.max)
        try:
            VarDbl(maxV, maxU) + VarDbl(maxV, maxU)
            self.fail("value overflow")
        except ValueException:
            pass
        except BaseException as ex:
            self.fail(ex)

        try:
            VarDbl(maxV, maxU) - VarDbl(maxV, maxU)
            self.fail("variance overflow")
        except UncertaintyException:
            pass
        except BaseException as ex:
            self.fail(ex)


class TestMultiply (unittest.TestCase):

    def testZero(self):
        validate(self, VarDbl(0) * VarDbl(1), 0, 0)
        validate(self, VarDbl(0) * 1, 0, 0)
        validate(self,  1 * VarDbl(0), 0, 0)

        validate(self, VarDbl(0) * VarDbl(1.0), 0, 0)
        validate(self, VarDbl(0) * 1.0, 0, 0)
        validate(self, 1.0 * VarDbl(0), 0, 0)

        validate(self, VarDbl(0, 1e-3) * VarDbl(2.0), 0, 2e-3)
        validate(self, VarDbl(0, 1e-3) * 2.0, 0, 2e-3)
        validate(self, 2.0 * VarDbl(0, 1e-3), 0, 2e-3)

        validate(self, VarDbl(0) * VarDbl(2.0, 1e-3), 0, 0)
        validate(self, VarDbl(2, 1e-3) * 0.0, 0, 0)
        validate(self, 0.0 * VarDbl(2, 1e-3), 0, 0)

        validate(self, VarDbl(0, 1e-3) * VarDbl(2.0), 0, 2e-3)
        validate(self, VarDbl(0, 1e-3) * 2.0, 0, 2e-3)
        validate(self, 2.0 * VarDbl(0, 1e-3), 0, 2e-3)

        validate(self, VarDbl(0, 1e-3) * VarDbl(2.0, 1e-2), 0, math.sqrt(4 + 1e-4)*1e-3)

    def testOne(self):
        validate(self, VarDbl(-1) * VarDbl(2), -2, 0)
        validate(self, VarDbl(-1) * 2, -2, 0)
        validate(self, 2 * VarDbl(-1), -2, 0)

        uncertainty = VarDbl.ulp(2.0)
        validate(self, VarDbl(-1.0) * VarDbl(2), -2, uncertainty)
        validate(self, VarDbl(-1.0) * 2, -2, uncertainty)
        validate(self, 2 * VarDbl(-1.0), -2, uncertainty)

        validate(self, VarDbl(-1.0) * VarDbl(2.0), -2, math.sqrt(2)*uncertainty)
        validate(self, VarDbl(-1.0, 1e-3) * VarDbl(2.0), -2, 2e-3)
        validate(self, VarDbl(-1.0) * VarDbl(2.0, 1e-3), -2, 1e-3)
        validate(self, VarDbl(-1.0, 1e-3) * VarDbl(2.0, 1e-3), -2, math.sqrt(5 + 1e-6) * 1e-3)

    def testException(self):
        maxV = math.sqrt(sys.float_info.max) * 1.001
        maxU = math.sqrt(sys.float_info.max) * 0.5
        try:
            VarDbl(maxV, maxU) * VarDbl(maxV, maxU)
            self.fail("value overflow")
        except ValueException:
            pass
        except BaseException as ex:
            self.fail(ex)

        try:
            VarDbl(1, maxU) * VarDbl(1, maxU)
            self.fail("variance overflow")
        except UncertaintyException:
            pass
        except BaseException as ex:
            self.fail(ex)

    def testSinUncertainty(self):
        res = VarDbl(1.1984225887068295e-05, 1.9530084268958033e-12) **2
        res2 = VarDbl(9.9999999992818911e-01, 9.3621846317751177e-17) **2
        res += res2
        res -= 1
        self.assertAlmostEqual(res.uncertainty(), 2.0337047278745187e-16, delta=4e-16)
        self.assertAlmostEqual(res.value(), -1.1102230246251565e-16, delta=4e-16)

        res = VarDbl(1.3295362544208494e-01, 2.7534451894481781e-16) **2
        res2 = VarDbl(9.9112225960363021e-01, 9.7483093233450446e-17) **2
        res += res2
        res -= 1
        self.assertAlmostEqual(res.uncertainty(), 2.1636206136563459e-16, delta=4e-16)
        self.assertEqual(res.value(), 0)


class TestPower (unittest.TestCase):
    def test_0(self):
        res = VarDbl(0, 1/1024) ** 0
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

        res = VarDbl(1, 1/1024) ** 0
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

    def test_1(self):
        res = VarDbl(-1, 1/1024) ** 1
        self.assertEqual(res.value(), -1)
        self.assertEqual(res.uncertainty(), 1/1024)


    def test_2(self):
        DELTA_VALUE = 3e-7
        DELTA_UNCERTAINTY = 2e-6
        
        for i in range(-16, 0):
            dx = math.pow(10, i)
            res = VarDbl(0, dx) ** 2
            try:
                self.assertAlmostEqual(res.value(), dx**2, delta=DELTA_VALUE)
                self.assertAlmostEqual(res.uncertainty(), math.sqrt(2)*(dx**2), delta=DELTA_UNCERTAINTY*4)
            except AssertionError as ex:
                raise ex        
        try:
            res = VarDbl(0, 1/8) ** (2 - 1e-9)
        except ZeroDivisionError:
            pass

        res = VarDbl(-1, 1/8) ** 2
        self.assertAlmostEqual(res.value(), 1 + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4 + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)
        try:
            res = VarDbl(-1, 1/8) ** (2 - 1e-9)
        except ValueError:
            pass

        res = VarDbl(1, 1/8) ** 2
        self.assertAlmostEqual(res.value(), 1 + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(1**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(1, 1/8) ** (2 - 1e-9)
        self.assertAlmostEqual(res.value(), 1 + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(1**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(1, 1/8) ** (2 + 1e-9)
        self.assertAlmostEqual(res.value(), 1 + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(1**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)

        res = VarDbl(2, 1/8) ** 2
        self.assertAlmostEqual(res.value(), (2**2) + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(2**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4*2)
        res = VarDbl(2, 1/8) ** (2 - 1e-9)
        self.assertAlmostEqual(res.value(), (2**2) + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(2**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4*2)
        res = VarDbl(2, 1/8) ** (2 + 1e-9)
        self.assertAlmostEqual(res.value(), (2**2) + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(2**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4*2)

        res = VarDbl(0.5, 1/8) ** 2
        self.assertAlmostEqual(res.value(), (0.5**2) + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(0.5**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(0.5, 1/8) ** (2 - 1e-9)
        self.assertAlmostEqual(res.value(), (0.5**2) + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(0.5**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(0.5, 1/8) ** (2 + 1e-9)
        self.assertAlmostEqual(res.value(), (0.5**2) + 1/(8**2), delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(4*(0.5**2) + 2/(8**2))/8, delta=DELTA_UNCERTAINTY*4)

    def test_3(self):
        DELTA_VALUE = 7e-7
        DELTA_UNCERTAINTY = 5e-6

        res = VarDbl(0, 1/8) ** 3
        self.assertAlmostEqual(res.value(), 0)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(15)/8**3, delta=3e-6)
        try:
            res = VarDbl(0, 1/8) ** (2 - 1e-9)
        except ZeroDivisionError:
            pass

        res = VarDbl(-1, 1/8) ** 3
        self.assertAlmostEqual(res.value(), -1 - 3/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9 + 36/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        res = VarDbl(1, 1/8) ** 3
        self.assertAlmostEqual(res.value(), 1 + 3/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9 + 36/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        res = VarDbl(1, 1/8) ** (3 - 1e-9)
        self.assertAlmostEqual(res.value(), 1 + 3/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9 + 36/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        res = VarDbl(1, 1/8) ** (3 + 1e-9)
        self.assertAlmostEqual(res.value(), 1 + 3/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9 + 36/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        try:
            res = VarDbl(-1, 1/8) ** (3 - 1e-9)
        except ValueError:
            pass

        res = VarDbl(-2, 1/8) ** 3
        self.assertAlmostEqual(res.value(), (-2**3) - 3*2/8**2, delta=DELTA_VALUE*2)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(2**4) + 36*(2**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(2, 1/8) ** 3
        self.assertAlmostEqual(res.value(), (2**3) + 3*2/8**2, delta=DELTA_VALUE*2)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(2**4) + 36*(2**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(2, 1/8) ** (3 - 1e-9)
        self.assertAlmostEqual(res.value(), (2**3) + 3*2/8**2, delta=DELTA_VALUE*2)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(2**4) + 36*(2**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY*4)
        res = VarDbl(2, 1/8) ** (3 + 1e-9)
        self.assertAlmostEqual(res.value(), (2**3) + 3*2/8**2, delta=DELTA_VALUE*2)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(2**4) + 36*(2**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY*4)

        res = VarDbl(-0.5, 1/8) ** 3
        self.assertAlmostEqual(res.value(), (-0.5**3) - 3*0.5/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(0.5**4) + 36*(0.5**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        res = VarDbl(0.5, 1/8) ** 3
        self.assertAlmostEqual(res.value(), (0.5**3) + 3*0.5/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(0.5**4) + 36*(0.5**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        res = VarDbl(0.5, 1/8) ** (3 - 1e-9)
        self.assertAlmostEqual(res.value(), (0.5**3) + 3*0.5/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(0.5**4) + 36*(0.5**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)
        res = VarDbl(0.5, 1/8) ** (3 + 1e-9)
        self.assertAlmostEqual(res.value(), (0.5**3) + 3*0.5/8**2, delta=DELTA_VALUE)
        self.assertAlmostEqual(res.uncertainty(), math.sqrt(9*(0.5**4) + 36*(0.5**2)/8**2 + 15/8**4)/8, delta=DELTA_UNCERTAINTY)


@unittest.skip('Not test')
class TestConvergence (unittest.TestCase):
    EDGE_HEADER = 'Edge Value\tEdge Uncertainty\tValue\tUncertainty\tException\n'

    def test_power(self):
        with open('./Python/Output/PowerAtOneEdge.txt', 'w') as f:
            f.write(TestConvergence.EDGE_HEADER)
            for i in range(-60, 60, 2):
                exp = i/10
                excpt = None
                for j in range(300,100,-1):
                    dx = j/1000
                    try:
                        res = VarDbl(1, dx)**exp
                    except BaseException as ex:
                        excpt = ex
                        continue
                    break
                if excpt:
                    f.write(f'{exp}\t{dx}\t{res.value()}\t{res.uncertainty()}\t{excpt}\n')
                    f.flush()

    def test_sin(self):
        with open('./Python/Output/SinEdge.txt', 'w') as f:
            f.write(TestConvergence.EDGE_HEADER)
            for i in range(-16, 17):
                rad = i/16 * math.pi
                excpt = None
                for j in range(150,50,-1):
                    dx = j/100
                    try:
                        res = VarDbl.sin(VarDbl(rad, dx))
                    except BaseException as ex:
                        excpt = ex
                        continue
                    break
                if excpt:
                    f.write(f'{rad}\t{dx}\t{res.value()}\t{res.uncertainty()}\t{excpt}\n')
                    f.flush()

    def test_log(self):
        with open('./Python/Output/LogEdge.txt', 'w') as f:
            f.write(TestConvergence.EDGE_HEADER)
            for i in range(-15, 17):
                x = i/16 + 1
                excpt = None
                for j in range(250,0,-1):
                    dx = j/100
                    try:
                        res = VarDbl.log(VarDbl(x, dx))
                    except BaseException as ex:
                        excpt = ex
                        continue
                    break
                if excpt and j:
                    f.write(f'{x}\t{dx}\t{res.value()}\t{res.uncertainty()}\t{excpt}\n')
                    f.flush()


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
        EXP = 2
        with open('./Python/Output/SquareAtZero.txt', 'w') as f:
            TestConvergence.writePowerHeader(f)
            for i in range(-4, 3):
                TestConvergence.calc_power(0, math.pow(10, i), EXP, f)
            for i in (1, EXP, 5):
                TestConvergence.calc_power(0, i * 1e3, EXP, f)

    def test_invesion_at_one(self):
        EXP = -1
        with self.assertRaises(NotMonotonicException):
            VarDbl(1, 0.1999872)**EXP
        res = VarDbl(1, 0.1999871)**EXP
        self.assertAlmostEqual(res.value(), 1.0462353179001609)
        self.assertAlmostEqual(res.uncertainty(), 0.24996133298478024)

        with open('./Python/Output/InversionAtOne.txt', 'w') as f:
            TestConvergence.writePowerHeader(f)
            for j in (0.1999871, 0.1):
                TestConvergence.calc_power(1, j, EXP, f) 
            for i in range(-4, -1):
                pw = math.pow(10, i)
                for x in (1, 2, 5):
                    TestConvergence.calc_power(1, x*pw, EXP, f) 
            TestConvergence.calc_power(1, 0.01, EXP, f)

    def test_square_root_at_one(self):
        EXP = 0.5
        with self.assertRaises(NotMonotonicException):
            VarDbl(1, 0.213)**EXP
        res = VarDbl(1, 0.212)**EXP
        self.assertAlmostEqual(res.value(), 0.9941081575421067)
        self.assertAlmostEqual(res.uncertainty(), 0.1083926741188568)

        with open('./Python/Output/SquareRootAtOne.txt', 'w') as f:
            TestConvergence.writePowerHeader(f)
            for i in range(-4, -1):
                pw = math.pow(10, i)
                for x in (1, 2, 5):
                    TestConvergence.calc_power(1, x*pw, EXP, f) 
            for x in (0.212, 0.2, 0.1):
                TestConvergence.calc_power(1, x, EXP, f) 

    def test_cubic_square_root_at_one(self):
        EXP = 1.5
        with self.assertRaises(NotReliableException):
            VarDbl(1, 0.244)**EXP
        res = VarDbl(1, 0.243)**EXP    # random.gauss(1, 0.243) will generate negative value
        self.assertAlmostEqual(res.value(), 1.0768954055344468)
        self.assertAlmostEqual(res.uncertainty(), 0.14527106865991737)

    def test_Inversion_square_one(self):
        EXP = -2
        with self.assertRaises(NotMonotonicException):
            VarDbl(1, 0.195)**EXP
        res = VarDbl(1, 0.194)**EXP    # random.gauss(1, 0.243) will generate negative value
        self.assertAlmostEqual(res.value(), 1.1436033475936196)
        self.assertAlmostEqual(res.uncertainty(), 0.7181123505697333)


    @staticmethod
    def writeHeader(f, divids:int=5, devs:int=3):
        f.write('Input Value\tInput Uncertainty'
                '\tValue\tUncertainty\tMean\tDeviation'
                '\tNormalized Error Mean\tNormalized Error Deviation\tLess\tMore')
        histo = Histo(divids, devs)
        for bucket in histo.buckets():
            f.write(f'\t{bucket:.1f}')
        f.write('\n')

    @staticmethod
    def calc(varFunc, mathFunc, x:float, dx:float, f):
        if dx <= 0:
            return
        CNT = 10000
        try:
            var = varFunc(VarDbl(x, dx))
        except BaseException as ex:
            raise ex
        varUnc = var.uncertainty()
        stat = Stat()
        histo = Histo(5, 3)
        for j in range(CNT):
            val = mathFunc(random.gauss(x, dx))
            try:
                stat.accum( val )
            except BaseException as ex:
                raise ex
            histo.accum( (val - var.value())/varUnc )
        f.write(f'{x}\t{dx}\t{var.value()}\t{var.uncertainty()}'
                f'\t{stat.mean()}\t{stat.dev()}\t{histo.stat().mean()}\t{histo.stat().dev()}\t{histo.less()}\t{histo.more()}')
        cnt = stat.count() - histo.less() - histo.more()
        for h in histo.histogram():
            f.write(f'\t{h/cnt}')
        f.write('\n')

    def test_sin_pi(self):
        X = 0
        with self.assertRaises(NotMonotonicException):
            VarDbl.sin(VarDbl(X, 0.734))
        res = VarDbl.sin(VarDbl(X, 0.733))    # random.gauss(1, 0.243) will generate negative value
        self.assertAlmostEqual(res.value(), 0)
        self.assertAlmostEqual(res.uncertainty(), 0.573828274225596)

        with open('./Python/Output/SinAtPI.txt', 'w') as f:
            TestConvergence.writeHeader(f)
            for i in range(-4, -1):
                pw = math.pow(10, i)
                for x in (1, 2, 5):
                    TestConvergence.calc(VarDbl.sin, math.sin, X, x*pw, f) 
            for x in (0.733, 0.5, 0.2, 0.1):
                TestConvergence.calc(VarDbl.sin, math.sin, X, x, f) 

    def test_sin_half_pi(self):
        X = math.pi/2
        with self.assertRaises(NotMonotonicException):
            VarDbl.sin(VarDbl(X, 0.574))
        res = VarDbl.sin(VarDbl(X, 0.573))
        self.assertAlmostEqual(res.value(), 0.848603461805469)
        self.assertAlmostEqual(res.uncertainty(), 0.19789632177500516)

        with open('./Python/Output/SinAtHalfPI.txt', 'w') as f:
            TestConvergence.writeHeader(f)
            for i in range(-4, -1):
                pw = math.pow(10, i)
                for x in (1, 2, 5):
                    TestConvergence.calc(VarDbl.sin, math.sin, X, x*pw, f) 
            for x in (0.573, 0.5, 0.2, 0.1):
                TestConvergence.calc(VarDbl.sin, math.sin, X, x, f) 

    def test_sin_quarter_pi(self):
        X = math.pi/4
        with self.assertRaises(NotMonotonicException):
            VarDbl.sin(VarDbl(X, 1.025))
        res = VarDbl.sin(VarDbl(X, 1.024))
        self.assertAlmostEqual(res.value(), 0.41859090238967567)
        self.assertAlmostEqual(res.uncertainty(), 0.5698961804018301)

        with open('./Python/Output/SinAtQuarterPI.txt', 'w') as f:
            TestConvergence.writeHeader(f)
            for i in range(-4, 0):
                pw = math.pow(10, i)
                for x in (1, 2, 5):
                    TestConvergence.calc(VarDbl.sin, math.sin, X, x*pw, f) 
            for x in (1.024, 1):
                TestConvergence.calc(VarDbl.sin, math.sin, X, x, f) 





class TestDivideBy (unittest.TestCase):

    def testVarDblByVarDblOne(self):
        validate(self, VarDbl(0) / VarDbl(1), 0, 0)
        validate(self, VarDbl(1) / VarDbl(1), 1, VarDbl.ulp(1))
        validate(self, VarDbl(-1) / VarDbl(1), -1, VarDbl.ulp(1))
        validate(self, VarDbl(1) / VarDbl(-1), -1, VarDbl.ulp(1))
        validate(self, VarDbl(-1) / VarDbl(-1), 1, VarDbl.ulp(1))
        validate(self, VarDbl(2) / VarDbl(1), 2, VarDbl.ulp(2))

        validate(self, VarDbl(1.0) / VarDbl(1), 1, math.sqrt(2)*VarDbl.ulp(1))
        validate(self, VarDbl(2.0) / VarDbl(1), 2, math.sqrt(2)*VarDbl.ulp(2))
        validate(self, VarDbl(0.5) / VarDbl(1), 0.5, math.sqrt(2)*VarDbl.ulp(0.5))

        validate(self, VarDbl(1.0) / VarDbl(1.0), 1, math.sqrt(3)*VarDbl.ulp(1), deltaUncertainty=5.30e-22)
        validate(self, VarDbl(2.0) / VarDbl(1.0), 2, math.sqrt(3)*VarDbl.ulp(2), deltaUncertainty=1.06e-21)
        validate(self, VarDbl(0.5) / VarDbl(1.0), 0.5, math.sqrt(3)*VarDbl.ulp(0.5), deltaUncertainty=2.65e-21)

        validate(self, VarDbl(0, 1e-3) / VarDbl(1), 0, 1e-3)
        validate(self, VarDbl(1, 1e-3) / VarDbl(1), 1, 1e-3)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(1), -1, 1e-3)
        validate(self, VarDbl(1, 1e-3) / VarDbl(-1), -1, 1e-3)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(-1), 1, 1e-3)
        validate(self, VarDbl(2, 1e-3) / VarDbl(1), 2, 1e-3)
        validate(self, VarDbl(0.5, 1e-3) / VarDbl(1), 0.5, 1e-3)

        validate(self, VarDbl(0) / VarDbl(1, 1e-3), 0, 0)
        validate(self, VarDbl(1) / VarDbl(1, 1e-3), 1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(-1) / VarDbl(1, 1e-3), -1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(1) / VarDbl(-1, 1e-3), -1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(-1) / VarDbl(-1, 1e-3), 1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(2) / VarDbl(1, 1e-3), 2, 0.002, deltaValue=2e-6, deltaUncertainty=6.4e-9)
        validate(self, VarDbl(0.5) / VarDbl(1, 1e-3), 0.5, 0.0005, deltaValue=5e-7, deltaUncertainty=4.0e-9)

        validate(self, VarDbl(0, 1e-3) / VarDbl(1, 1e-3), 0, 1e-3, deltaUncertainty=1.5e-9)
        validate(self, VarDbl(1, 1e-3) / VarDbl(1, 1e-3), 1, 1e-3*math.sqrt(2), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(1, 1e-3), -1, 1e-3*math.sqrt(2), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(1, 1e-3) / VarDbl(-1, 1e-3), -1, 1e-3*math.sqrt(2), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(-1, 1e-3), 1, 1e-3*math.sqrt(2), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(2, 1e-3) / VarDbl(1, 1e-3), 2, 1e-3*math.sqrt(5), deltaValue=2e-6, deltaUncertainty=5.0e-9)
        validate(self, VarDbl(0.5, 1e-3) / VarDbl(1, 1e-3), 0.5, 0.0005*math.sqrt(5), deltaValue=5e-7, deltaUncertainty=1.2e-9)

    def testVarDblByFloatOne(self):
        validate(self, VarDbl(0)/ 1.0, 0, 0)
        validate(self, VarDbl(1)/ 1.0,   1, math.sqrt(2)*VarDbl.ulp(1), deltaUncertainty=6.5e-22)
        validate(self, VarDbl(-1)/ 1.0, -1, math.sqrt(2)*VarDbl.ulp(1), deltaUncertainty=6.5e-22)
        validate(self, VarDbl(1)/ -1.0, -1, math.sqrt(2)*VarDbl.ulp(1), deltaUncertainty=6.5e-22)
        validate(self, VarDbl(-1)/ -1.0, 1, math.sqrt(2)*VarDbl.ulp(1), deltaUncertainty=6.5e-22)
        validate(self, VarDbl(2)/ 1.0, 2,   2*math.sqrt(2)*VarDbl.ulp(1), deltaUncertainty=1.3e-21)

        validate(self, VarDbl(1.0)/ 1.0, 1,   math.sqrt(3)*VarDbl.ulp(1), deltaUncertainty=1.3e-21)
        validate(self, VarDbl(2.0)/ 1.0, 2,   math.sqrt(3)*VarDbl.ulp(2), deltaUncertainty=1.3e-21)
        validate(self, VarDbl(0.5)/ 1.0, 0.5, math.sqrt(3)*VarDbl.ulp(0.5), deltaUncertainty=2.7e-22)

        validate(self, VarDbl(0, 1e-3)/ 1.0, 0, 1e-3)
        validate(self, VarDbl(1, 1e-3)/ 1.0, 1, 1e-3)
        validate(self, VarDbl(-1, 1e-3)/ 1.0, -1, 1e-3)
        validate(self, VarDbl(1, 1e-3)/ -1.0, -1, 1e-3)
        validate(self, VarDbl(-1, 1e-3)/ -1.0, 1, 1e-3)
        validate(self, VarDbl(2, 1e-3)/ 1.0, 2, 1e-3)
        validate(self, VarDbl(0.5, 1e-3)/ 1.0, 0.5, 1e-3)

    def testFloatByVarDblOne(self):
        validate(self, 0 / VarDbl(1), 0, 0)
        validate(self, 1 / VarDbl(1), 1, VarDbl.ulp(1))
        validate(self, -1 / VarDbl(1), -1, VarDbl.ulp(1))
        validate(self, 1 / VarDbl(-1), -1, VarDbl.ulp(1))
        validate(self, -1 / VarDbl(-1), 1, VarDbl.ulp(1))
        validate(self, 2 / VarDbl(1), 2, VarDbl.ulp(2))
        validate(self, 0.5 / VarDbl(1), 0.5, math.ulp(0.5)/math.sqrt(3/2))

        validate(self, 0 / VarDbl(1, 1e-3), 0, 0)
        validate(self, 1 / VarDbl(1, 1e-3), 1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, -1 / VarDbl(1, 1e-3), -1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, 1 / VarDbl(-1, 1e-3), -1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, -1 / VarDbl(-1, 1e-3), 1, 1e-3, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, 2 / VarDbl(1, 1e-3), 2, 0.002, deltaValue=2e-6, deltaUncertainty=6.4e-9)
        validate(self, 0.5 / VarDbl(1, 1e-3), 0.5, 0.0005, deltaValue=5e-7, deltaUncertainty=4.0e-9)

    def testVarDblByVarDblZero(self):
        try:
            VarDbl(0) / VarDbl(0)
            self.fail()
        except ValueError:
            pass
        except Exception as ex:
            self.fail(ex)

        try:
            VarDbl(0) / VarDbl(0, 1)
            self.fail()
        except ValueError:
            pass
        except Exception as ex:
            self.fail(ex)

        try:
            VarDbl(0) / VarDbl(0.1, 1)
            self.fail()
        except NotReliableException:
            pass
        except NotMonotonicException:
            pass
        except Exception as ex:
            self.fail(ex)

    def testVarDblByVarDblTwo(self):
        validate(self, VarDbl(0) / VarDbl(2), 0, 0)
        validate(self, VarDbl(1) / VarDbl(2), 0.5, VarDbl.ulp(0.5))
        validate(self, VarDbl(-1) / VarDbl(2), -0.5, VarDbl.ulp(0.5))
        validate(self, VarDbl(1) / VarDbl(-2), -0.5, VarDbl.ulp(0.5))
        validate(self, VarDbl(-1) / VarDbl(-2), 0.5, VarDbl.ulp(0.5))
        validate(self, VarDbl(2) / VarDbl(2), 1, VarDbl.ulp(1))
        validate(self, VarDbl(0.5) / VarDbl(2), 0.25, math.sqrt(2)*VarDbl.ulp(0.25))

        validate(self, VarDbl(0, 1e-3) / VarDbl(2), 0, 5e-4)
        validate(self, VarDbl(1, 1e-3) / VarDbl(2), 0.5, 5e-4)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(2), -0.5, 5e-4)
        validate(self, VarDbl(1, 1e-3) / VarDbl(-2), -0.5, 5e-4)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(-2), 0.5, 5e-4)
        validate(self, VarDbl(2, 1e-3) / VarDbl(2), 1, 5e-4)
        validate(self, VarDbl(0.5, 1e-3) / VarDbl(2), 0.25, 5e-4)

        validate(self, VarDbl(0) / VarDbl(2, 1e-3), 0, 0)
        validate(self, VarDbl(1) / VarDbl(2, 1e-3), 0.5, 2.5e-4, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(-1) / VarDbl(2, 1e-3), -0.5, 2.5e-4, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(1) / VarDbl(-2, 1e-3), -0.5, 2.5e-4, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(-1) / VarDbl(-2, 1e-3), 0.5, 2.5e-4, deltaValue=1e-6, deltaUncertainty=4.0e-9)
        validate(self, VarDbl(2) / VarDbl(2, 1e-3), 1, 5e-4, deltaValue=2e-6, deltaUncertainty=6.4e-9)
        validate(self, VarDbl(0.5) / VarDbl(2, 1e-3), 0.25, 1.25e-4, deltaValue=5e-7, deltaUncertainty=4.0e-9)

        validate(self, VarDbl(0, 1e-3) / VarDbl(2, 1e-3), 0, 5e-4, deltaUncertainty=1.5e-9)
        validate(self, VarDbl(1, 1e-3) / VarDbl(2, 1e-3), 0.5, 2.5e-4*math.sqrt(5), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(2, 1e-3), -0.5, 2.5e-4*math.sqrt(5), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(1, 1e-3) / VarDbl(-2, 1e-3), -0.5, 2.5e-4*math.sqrt(5), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(-1, 1e-3) / VarDbl(-2, 1e-3), 0.5, 2.5e-4*math.sqrt(5), deltaValue=1e-6, deltaUncertainty=1.2e-9)
        validate(self, VarDbl(2, 1e-3) / VarDbl(2, 1e-3), 1, 1e-3*math.sqrt(0.5), deltaValue=2e-6, deltaUncertainty=5.0e-9)
        validate(self, VarDbl(0.5, 1e-3) / VarDbl(2, 1e-3), 0.25, 2.5e-4*math.sqrt(4.25), deltaValue=5e-7, deltaUncertainty=1.2e-9)





class TestCompare (unittest.TestCase):

    def testErf(self):
        self.assertAlmostEqual(0, math.erf(0), delta=1e-3)
        self.assertAlmostEqual(0.8427, math.erf(1), delta=1e-3)
        self.assertAlmostEqual(-0.8427, math.erf(-1), delta=1e-3)
        # Normal bounds
        self.assertAlmostEqual(0.6827, math.erf(1/math.sqrt(2)), delta=1e-3)
        self.assertAlmostEqual(0.9545, math.erf(2/math.sqrt(2)), delta=1e-3)
        self.assertAlmostEqual(0.9973, math.erf(3/math.sqrt(2)), delta=1e-3)
        
    def testBindingForEqual(self):
        self.assertAlmostEqual(0.5, math.erf(VarDbl.BINDING_FOR_EQUAL/math.sqrt(2)), delta=1e-7)
        self.assertAlmostEqual(VarDbl.BINDING_FOR_EQUAL, math.erfc(0.5)*math.sqrt(2), delta=5e-3)

    def testVarDbl(self):
        v1 = VarDbl(1.000, 0.002)
        v2 = VarDbl(1.001, 1e-3)
        self.assertTrue(v1 == v2)
        self.assertFalse(v1 != v2)
        self.assertFalse(v1 < v2)
        self.assertTrue(v1 <= v2)
        self.assertFalse(v1 > v2)
        self.assertTrue(v1 >= v2)

        v3 = VarDbl(1.002, 1e-3)
        self.assertFalse(v1 == v3)
        self.assertTrue(v1 != v3)
        self.assertTrue(v1 < v3)
        self.assertTrue(v1 <= v3)
        self.assertFalse(v1 > v3)
        self.assertFalse(v1 >= v3)
        self.assertFalse(v3 < v1)
        self.assertFalse(v3 <= v1)
        self.assertTrue(v3 > v1)
        self.assertTrue(v3 >= v1)

    def testFloat(self):
        v1 = VarDbl(1.000, 0.002)
        v2 = 1.001
        self.assertTrue(v1 == v2)
        self.assertFalse(v1 != v2)
        self.assertFalse(v1 < v2)
        self.assertTrue(v1 <= v2)
        self.assertFalse(v1 > v2)
        self.assertTrue(v1 >= v2)

        self.assertTrue(v2 == v1)
        self.assertFalse(v2 != v1)
        self.assertFalse(v2 < v1)
        self.assertTrue(v2 <= v1)
        self.assertFalse(v2 > v1)
        self.assertTrue(v2 >= v1)

        v1 = VarDbl(1.000, 0.002)
        v3 = 1.002
        self.assertFalse(v1 == v3)
        self.assertTrue(v1 != v3)
        self.assertTrue(v1 < v3)
        self.assertTrue(v1 <= v3)
        self.assertFalse(v1 > v3)
        self.assertFalse(v1 >= v3)
        self.assertFalse(v3 < v1)
        self.assertFalse(v3 <= v1)
        self.assertTrue(v3 > v1)
        self.assertTrue(v3 >= v1)

 






if __name__ == '__main__':
    unittest.main()