from contextlib import AbstractContextManager
import math
import pickle
import unittest
import sys

from varDbl import VarDbl, InitException, InitException, validate
from taylor import NotReliableException, NotMonotonicException


class TestInit (unittest.TestCase):

    def testInt(self):
        validate(self, VarDbl(0), 0, 0)
        validate(self, VarDbl(1), 1, 0)
        validate(self, VarDbl(-1), -1, 0)


    def testLargeInt(self):
        validate(self, VarDbl(VarDbl.DOUBLE_MAX_SIGNIFICAND), VarDbl.DOUBLE_MAX_SIGNIFICAND, 0)
        validate(self, VarDbl(-VarDbl.DOUBLE_MAX_SIGNIFICAND), -VarDbl.DOUBLE_MAX_SIGNIFICAND, 0)

        # lost resolution
        f = float(VarDbl.DOUBLE_MAX_SIGNIFICAND + 1)
        validate(self, VarDbl(VarDbl.DOUBLE_MAX_SIGNIFICAND + 1), f, 0.5)
        # no lost of resolution
        f = float(1 << 54)
        validate(self, VarDbl(1 << 54), f, 0)

    def testInit(self):
        validate(self, VarDbl(-1, 0), -1, 0)
        validate(self, VarDbl(-1, 1), -1, 1)
        validate(self, VarDbl(-1, -1), -1, 1)

    def testValueException(self):
        with self.assertRaises(InitException):
            VarDbl(float('nan'))
        with self.assertRaises(InitException):
            VarDbl(float('inf'))

    def testVarianceException(self):
        with self.assertRaises(InitException):
            VarDbl(0, float('nan'))
        with self.assertRaises(InitException):
            VarDbl(0, float('inf'))
            
    def testUncertaintyRange(self):
        maxU = math.sqrt(sys.float_info.max)
        with self.assertRaises(InitException):
            VarDbl(0, maxU + VarDbl.ulp(maxU))
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
        validate(self, v2 + 2.0, 3, 0) 
        validate(self, 2.0 + v2, 3, 0) 
        v2 += 2.0
        validate(self, v2, 3, 0) 

    def testSubFloat(self):
        v1 = VarDbl(1, math.sqrt(2))
        validate(self, v1 - 2.0, -1, math.sqrt(2)) 
        v = 2.0 - v1
        validate(self, 2.0 - v1, 1, math.sqrt(2)) 
        v1 -= 2.0
        validate(self, v1, -1, math.sqrt(2)) 

        v2 = VarDbl(1.0)
        validate(self, v2 - 2.0, -1, 0) 
        validate(self, 2.0 - v2, 1, 0) 
        v2 -= 2.0
        validate(self, v2, -1, 0) 

    def testAddSubException(self):
        maxV = sys.float_info.max
        maxU = math.sqrt(sys.float_info.max)
        try:
            VarDbl(maxV, maxU) + VarDbl(maxV, maxU)
            self.fail("value overflow")
        except InitException:
            pass
        except BaseException as ex:
            self.fail(ex)

        try:
            VarDbl(maxV, maxU) - VarDbl(maxV, maxU)
            self.fail("variance overflow")
        except InitException:
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

        validate(self, VarDbl(-1.0) * VarDbl(2), -2, 0)
        validate(self, VarDbl(-1.0) * 2, -2, 0)
        validate(self, 2 * VarDbl(-1.0), -2, 0)

        validate(self, VarDbl(-1.0) * VarDbl(2.0), -2, 0)
        validate(self, VarDbl(-1.0, 1e-3) * VarDbl(2.0), -2, 2e-3)
        validate(self, VarDbl(-1.0) * VarDbl(2.0, 1e-3), -2, 1e-3)
        validate(self, VarDbl(-1.0, 1e-3) * VarDbl(2.0, 1e-3), -2, math.sqrt(5 + 1e-6) * 1e-3)

    def testTwoWithHalf(self):
        validate(self, VarDbl(0.5) * VarDbl(2), 1, 0)
        validate(self, VarDbl(0.5, 0) * 2, 1,  0)
        validate(self, 0.5 * VarDbl(-2), -1, 0)

        validate(self, VarDbl(-0.5, 1e-3) * VarDbl(2.0), -1, 2e-3)
        validate(self, VarDbl(-0.5) * VarDbl(2.0, 1e-3), -1, 0.5e-3)
        validate(self, VarDbl(-0.5, 1e-3) * VarDbl(2.0, 1e-3), -1, math.sqrt(4.25 + 1e-6) * 1e-3)

    def testException(self):
        maxV = math.sqrt(sys.float_info.max) * 1.001
        maxU = math.sqrt(sys.float_info.max) * 0.5
        try:
            VarDbl(maxV, maxU) * VarDbl(maxV, maxU)
            self.fail("value overflow")
        except InitException:
            pass
        except BaseException as ex:
            self.fail(ex)

        try:
            VarDbl(1, maxU) * VarDbl(1, maxU)
            self.fail("variance overflow")
        except InitException:
            pass
        except BaseException as ex:
            self.fail(ex)

    def testLargeDiff(self):
        self.assertEqual(13316075197586562, 64919121*205117922)
        self.assertEqual(13316075197586561, 159018721*83739041)
        self.assertEqual(1, 64919121*205117922 - 159018721*83739041)

        re12 = VarDbl(64919121) * VarDbl(205117922)
        self.assertEqual(13316075197586562.0, re12.value(), math.ulp(re12.value()))
        self.assertEqual(0, re12.uncertainty(), math.ulp(re12.uncertainty()))

        re34 = VarDbl(-159018721) * VarDbl(83739041)
        self.assertEqual( -13316075197586560.0, re34.value(), math.ulp(re34.value()))
        self.assertEqual(0.5, re34.uncertainty(), math.ulp(re34.uncertainty()))

        re = re12 + re34
        self.assertEqual(2, re.value(), math.ulp(re.value()))
        self.assertEqual(VarDbl.ulp(re12.value()), VarDbl.ulp(re34.value()), math.ulp(re.uncertainty()))
        self.assertEqual(0.5, re.uncertainty(), math.ulp(re.uncertainty()))

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




class TestDivideBy (unittest.TestCase):

    def testVarDblByVarDblOne(self):
        validate(self, VarDbl(0) / VarDbl(1), 0, 0)
        validate(self, VarDbl(1) / VarDbl(1), 1, 0)
        validate(self, VarDbl(-1) / VarDbl(1), -1, 0)
        validate(self, VarDbl(1) / VarDbl(-1), -1, 0)
        validate(self, VarDbl(-1) / VarDbl(-1), 1, 0)
        validate(self, VarDbl(2) / VarDbl(1), 2, 0)

        validate(self, VarDbl(1.0) / VarDbl(1), 1, 0)
        validate(self, VarDbl(2.0) / VarDbl(1), 2, 0)
        validate(self, VarDbl(0.5) / VarDbl(1), 0.5, 0)

        validate(self, VarDbl(1.0) / VarDbl(1.0), 1, 0)
        validate(self, VarDbl(2.0) / VarDbl(1.0), 2, 0)
        validate(self, VarDbl(0.5) / VarDbl(1.0), 0.5, 0)

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
        validate(self, VarDbl(1)/ 1.0,   1, 0)
        validate(self, VarDbl(-1)/ 1.0, -1, 0)
        validate(self, VarDbl(1)/ -1.0, -1, 0)
        validate(self, VarDbl(-1)/ -1.0, 1, 0)
        validate(self, VarDbl(2)/ 1.0, 2,   0)

        validate(self, VarDbl(1.0)/ 1.0, 1,   0)
        validate(self, VarDbl(2.0)/ 1.0, 2,   0)
        validate(self, VarDbl(0.5)/ 1.0, 0.5, 0)

        validate(self, VarDbl(0, 1e-3)/ 1.0, 0, 1e-3)
        validate(self, VarDbl(1, 1e-3)/ 1.0, 1, 1e-3)
        validate(self, VarDbl(-1, 1e-3)/ 1.0, -1, 1e-3)
        validate(self, VarDbl(1, 1e-3)/ -1.0, -1, 1e-3)
        validate(self, VarDbl(-1, 1e-3)/ -1.0, 1, 1e-3)
        validate(self, VarDbl(2, 1e-3)/ 1.0, 2, 1e-3)
        validate(self, VarDbl(0.5, 1e-3)/ 1.0, 0.5, 1e-3)

    def testFloatByVarDblOne(self):
        validate(self, 0 / VarDbl(1), 0, 0)
        validate(self, 1 / VarDbl(1), 1, 0)
        validate(self, -1 / VarDbl(1), -1, 0)
        validate(self, 1 / VarDbl(-1), -1, 0)
        validate(self, -1 / VarDbl(-1), 1, 0)
        validate(self, 2 / VarDbl(1), 2, 0)
        validate(self, 0.5 / VarDbl(1), 0.5, 0)

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
        validate(self, VarDbl(1) / VarDbl(2), 0.5, 0)
        validate(self, VarDbl(-1) / VarDbl(2), -0.5, 0)
        validate(self, VarDbl(1) / VarDbl(-2), -0.5, 0)
        validate(self, VarDbl(-1) / VarDbl(-2), 0.5, 0)
        validate(self, VarDbl(2) / VarDbl(2), 1, 0)
        validate(self, VarDbl(0.5) / VarDbl(2), 0.25, 0)

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