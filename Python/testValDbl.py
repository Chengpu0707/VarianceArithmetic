import math
import pickle
import unittest
import sys

from VarDbl import VarDbl, ValDblValueError, ValDblUncertaintyError

class TestInit (unittest.TestCase):

    def testInt(self):
        v = VarDbl(0)
        self.assertEqual(0, v.value())
        self.assertEqual(0, v.uncertainty())

        v = VarDbl(1)
        self.assertEqual(1, v.value())
        self.assertEqual(0, v.uncertainty())

        v = VarDbl(-1)
        self.assertEqual(-1, v.value())
        self.assertEqual(0, v.uncertainty())


    def testLargeInt(self):
        v = VarDbl((1 << 53))
        self.assertEqual((1 << 53), v.value())
        self.assertEqual(0, v.uncertainty())
        v = VarDbl(-(1 << 53))
        self.assertEqual(-(1 << 53), v.value())
        self.assertEqual(0, v.uncertainty())

        # lost resolution
        v = VarDbl((1 << 53) + 1)
        f = float((1 << 53) + 1)
        self.assertEqual(f, v.value())
        self.assertEqual(math.ulp(f), v.uncertainty())

    def testInit(self):
        v = VarDbl(-1, 0)
        self.assertEqual(-1, v.value())
        self.assertEqual(0, v.uncertainty())

        v = VarDbl(-1, 1)
        self.assertEqual(-1, v.value())
        self.assertEqual(1, v.uncertainty())

        v = VarDbl(-1, -1)
        self.assertEqual(-1, v.value())
        self.assertEqual(1, v.uncertainty())

    def failValDblValueError(self, value):
        try:
            VarDbl(value)
            self.fail(f'Init ValDbl with {value}')
        except ValDblValueError as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)

    def testValDblValueError(self):
        self.failValDblValueError(float('nan'))
        self.failValDblValueError(float('inf'))

    def failValDblUncertaintyError(self, value, uncertainty):
        try:
            VarDbl(value, uncertainty)
            self.fail(f'Init ValDbl with {value}~{uncertainty}')
        except ValDblUncertaintyError as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)

    def testValDblUncertaintyError(self):
        self.failValDblUncertaintyError(0, float('nan'))
        self.failValDblUncertaintyError(0, float('inf'))
            
    def testUncertaintyRange(self):
        maxU = math.sqrt(sys.float_info.max)
        minU = math.sqrt(math.ulp(sys.float_info.min))

        self.failValDblUncertaintyError(0, maxU + math.ulp(maxU))

        v = VarDbl(sys.float_info.max, maxU - math.ulp(maxU))
        self.assertEqual(sys.float_info.max, v.value())
        self.assertTrue(maxU - math.ulp(maxU) <= v.uncertainty() <= maxU + math.ulp(maxU))
        
        v = VarDbl(sys.float_info.min, minU)
        self.assertEqual(sys.float_info.min, v.value())
        self.assertTrue(minU - math.ulp(minU) <= v.uncertainty() <= minU + math.ulp(minU))

        v = VarDbl(sys.float_info.min, minU * 0.5)
        self.assertEqual(sys.float_info.min, v.value())
        self.assertAlmostEqual(0, v.uncertainty(), math.ulp(minU))


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
        self.assertEqual(v._value, vr._value, math.ulp(v._value))
        self.assertEqual(v._variance, vr._variance, math.ulp(v._value))




if __name__ == '__main__':
    unittest.main()