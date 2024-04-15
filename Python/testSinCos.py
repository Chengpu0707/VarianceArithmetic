'''
The resolution of floating input is a major source of error.
Folding input to [-pi/4, +pi/4] will not reduce the numerical error

The input should be index frequence to reduce the error.
'''


import math
import unittest

def sin(rad:float) -> float:
    quard = int(math.floor(rad * 4 / math.pi))
    if abs(quard % 2) == 1:
        quard += 1
    rem = rad - quard * math.pi / 4
    quard //= 2
    match quard % 4:
        case 0:
            return math.sin(rem)
        case 1:
            return math.cos(rem)
        case 2:
            return math.sin(-rem)
        case 3:
            return -math.cos(rem)
    
    

def cos(rad:float) -> float:
    rad = abs(rad)
    if rad <= math.pi / 4:
        return math.cos(rad)
    return sin(rad + math.pi/2)



class TestSinCos (unittest.TestCase):

    def validate(self, expected, result, delta=None):
        if delta is None:
            delta=math.ulp(expected)
        self.assertAlmostEqual(expected, result, delta=delta)

    def testRem(self):
        self.assertEqual(2, 10 % 8)
        self.assertEqual(6, -10 % 8)
        self.assertEqual(-2, -10 % -8)

    def testInQuard(self):
        self.validate(0, sin(0))
        self.validate( math.sin(3/16*math.pi), sin( 3/16*math.pi))
        self.validate(-math.sin(3/16*math.pi), sin(-3/16*math.pi)) 

        self.validate(1, cos(0))
        self.validate( math.cos(3/16*math.pi), cos( 3/16*math.pi))
        self.validate( math.cos(3/16*math.pi), cos(-3/16*math.pi)) 

    def testBoundary(self):
        self.assertEqual(math.sin(1/4*math.pi), math.cos(1/4*math.pi))
        self.assertEqual(math.sin(1/4*math.pi), math.cos(-1/4*math.pi))

        self.validate( math.sin(1/4*math.pi), sin( 1/4*math.pi))
        self.validate(-math.sin(1/4*math.pi), sin(-1/4*math.pi))

        self.validate(math.cos(1/4*math.pi), cos( 1/4*math.pi))
        self.validate(math.cos(1/4*math.pi), cos(-1/4*math.pi))

    def testInPeriod(self):
        self.validate( math.cos(3/16*math.pi), sin( 5/16*math.pi))
        self.validate(-math.cos(3/16*math.pi), sin(-5/16*math.pi))
        self.validate( math.cos(1/16*math.pi), sin( 7/16*math.pi))
        self.validate(-math.cos(1/16*math.pi), sin(-7/16*math.pi))
        self.validate( math.cos(1/16*math.pi), sin( 9/16*math.pi))
        self.validate(-math.cos(1/16*math.pi), sin(-9/16*math.pi))
        self.validate( math.cos(3/16*math.pi), sin( 11/16*math.pi))
        self.validate(-math.cos(3/16*math.pi), sin(-11/16*math.pi))

        self.validate( math.sin(3/16*math.pi), sin( 13/16*math.pi), 2.220446049250313e-16)
        self.validate(-math.sin(3/16*math.pi), sin(-13/16*math.pi), 2.220446049250313e-16)
        self.validate( math.sin(1/16*math.pi), sin( 15/16*math.pi), 2.220446049250313e-16)
        self.validate(-math.sin(1/16*math.pi), sin(-15/16*math.pi), 2.220446049250313e-16)
        self.validate(-math.sin(1/16*math.pi), sin( 17/16*math.pi), 2.220446049250313e-16)
        self.validate( math.sin(1/16*math.pi), sin(-17/16*math.pi), 2.220446049250313e-16)

    def testCompare(self):
        self.validate( math.sin(3/16*math.pi), math.sin( 13/16*math.pi))
        self.validate(-math.sin(3/16*math.pi), math.sin(-13/16*math.pi))
        self.validate( math.sin(1/16*math.pi), math.sin( 15/16*math.pi), 3.608224830031759e-16)
        self.validate(-math.sin(1/16*math.pi), math.sin(-15/16*math.pi), 3.608224830031759e-16)
        self.validate(-math.sin(1/16*math.pi), math.sin( 17/16*math.pi), 1.1102230246251565e-16)
        self.validate( math.sin(1/16*math.pi), math.sin(-17/16*math.pi), 1.1102230246251565e-16)


    def testError(self):
        self.validate(0, math.sin(0))
        self.validate(0, math.sin(   math.pi), 1.2246467991473532e-16)
        self.validate(0, math.sin(  -math.pi), 1.2246467991473532e-16)
        self.validate(0, math.sin( 2*math.pi), 2.4492935982947064e-16)
        self.validate(0, math.sin(-2*math.pi), 2.4492935982947064e-16)

        self.validate(0, math.cos(   math.pi/2),  6.123233995736766e-17)
        self.validate(0, math.cos(  -math.pi/2),  6.123233995736766e-17)
        self.validate(0, math.cos( 3*math.pi/2), 1.8369701987210297e-16)
        self.validate(0, math.cos(-3*math.pi/2), 1.8369701987210297e-16)


if __name__ == '__main__':
    unittest.main()