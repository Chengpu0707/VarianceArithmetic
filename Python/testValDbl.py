import math
import pickle
import unittest
import sys

from varDbl import VarDbl, ValueException, UncertaintyException

def validate(self, res, value, uncertainty):
    self.assertAlmostEqual(value, res.value(), delta=math.ulp(res.value()))
    self.assertAlmostEqual(uncertainty, res.uncertainty(), delta=math.ulp(res.uncertainty()))



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
        validate(self, VarDbl((1 << 53) + 1), f, VarDbl.DEVIATION_OF_LSB * math.ulp(f))
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
        self.failUncertaintyException(0, maxU + VarDbl.DEVIATION_OF_LSB * math.ulp(maxU))
        validate(self, VarDbl(0, maxU), 0, maxU) 
        
        minU = math.sqrt(VarDbl.DEVIATION_OF_LSB * math.ulp(sys.float_info.min))
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
        uncertainty = VarDbl.DEVIATION_OF_LSB * math.ulp(1) * math.sqrt(5)
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
        uncertainty = VarDbl.DEVIATION_OF_LSB * math.ulp(1) * math.sqrt(5)
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

        uncertainty = VarDbl.DEVIATION_OF_LSB * math.ulp(2.0)
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
        v2 = VarDbl(1.001, 0.001)
        self.assertTrue(v1 == v2)
        self.assertFalse(v1 != v2)
        self.assertFalse(v1 < v2)
        self.assertTrue(v1 <= v2)
        self.assertFalse(v1 > v2)
        self.assertTrue(v1 >= v2)

        v3 = VarDbl(1.002, 0.001)
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