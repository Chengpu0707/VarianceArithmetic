import math
import unittest

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

        v = VarDbl((1 << 53) + 1)
        f = float((1 << 53) + 1)
        self.assertEqual(f, v.value())
        self.assertEqual(math.ulp(f), v.uncertainty())


    def testValDblValueError(self):
        try:
            VarDbl(float('nan'))
            self.fail('Init ValDbl with nan')
        except ValDblValueError as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)
        try:
            VarDbl(float('inf'))
            self.fail('Init ValDbl with inf')
        except ValDblValueError as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)

    def testValDblUncertaintyError(self):
        try:
            VarDbl(0, float('nan'))
            self.fail('Init ValDbl with nan')
        except ValDblUncertaintyError as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)
        try:
            VarDbl(0, float('inf'))
            self.fail('Init ValDbl with inf')
        except ValDblUncertaintyError as ex:
            self.assertIsNotNone(ex.__traceback__)
        except BaseException as ex:
            self.fail(ex)


if __name__ == '__main__':
    unittest.main()