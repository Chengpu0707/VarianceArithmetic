from fractions import Fraction
import itertools
import numpy
import unittest

from matrix import permutSign, isSquareMatrix, createIntMatrix, addNoise, determinant
from varDbl import VarDbl


class TestPermumation (unittest.TestCase):
    def testPermut(self):
        '''
        The permutation is not [+, -]*
        '''
        self.assertTupleEqual(tuple(itertools.permutations(range(2), 2)), 
                              ((0,1), (1,0)))
        self.assertTupleEqual(tuple(itertools.permutations(range(3), 3)), 
                              ((0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)))
        self.assertTupleEqual(tuple(itertools.permutations(range(3), 2)), 
                              ((0,1), (0,2), (1,0), (1,2), (2,0), (2,1)))

    def testSign(self):
        self.assertEqual(permutSign((0,1)), 1)
        self.assertEqual(permutSign((1,0)), -1)

        self.assertEqual(permutSign((0,1,2)), 1)
        self.assertEqual(permutSign((0,2,1)), -1)
        self.assertEqual(permutSign((1,0,2)), -1)
        self.assertEqual(permutSign((1,2,0)), 1)
        self.assertEqual(permutSign((2,0,1)), 1)
        self.assertEqual(permutSign((2,1,0)), -1)


class TestSquareMatrix (unittest.TestCase):
    def testIsSquareMatrix(self):
        self.assertTrue( isSquareMatrix(tuple([tuple([VarDbl(1), 2.0]), tuple([3, VarDbl(4, 1)])])))

        self.assertFalse(isSquareMatrix(tuple([tuple([VarDbl(1), 2.0]), [3, VarDbl(4, 1)]])))
        self.assertFalse(isSquareMatrix(      [tuple([VarDbl(1), 2.0]), tuple([3, VarDbl(4, 1)])]))

        self.assertFalse(isSquareMatrix(tuple([tuple([VarDbl(1), 2.0]), tuple([3, None])])))
        self.assertFalse(isSquareMatrix(tuple([tuple([VarDbl(1), 2.0]), tuple([3, VarDbl(4, 1)])]), sType=(int, float)))

    def testCreateIntMatrix(self):
        '''
        Can not use numpy to hold int matrix and to calculate determinant
        '''
        matrix = createIntMatrix(2)
        self.assertTrue(isSquareMatrix(matrix, sType=(int,)))
        mat = numpy.matrix(matrix)
        self.assertEqual(mat.dtype, numpy.int32)
        self.assertEqual(type(numpy.linalg.det(mat)), numpy.float64)

    def testAddIntNoise(self):
        ssMatrix = addNoise(createIntMatrix(2), 1e-3)
        self.assertTrue(isSquareMatrix(ssMatrix, sType=(VarDbl,)))
        self.assertFalse([type(val.value()) for sMatrix in ssMatrix for val in sMatrix if type(val.value()) != int])

    def testAddFloatNoise(self):
        ssMatrix = addNoise(createIntMatrix(2), 1e-3, retainIntValue=False)
        self.assertTrue(isSquareMatrix(ssMatrix, sType=(VarDbl,)))
        self.assertFalse([type(val.value()) for sMatrix in ssMatrix for val in sMatrix if type(val.value()) != float])

class TestFractions (unittest.TestCase): 
    '''
    The promotion is int -> Fraction -> float
    '''
    def testFractions(self):
        self.assertEqual(Fraction(1,2) * Fraction(-4,3) + Fraction(5,6), Fraction(1,6))

    def testInt(self):
        self.assertEqual(Fraction(1,2) * (-3), Fraction(-3,2))

    def testInt(self):
        self.assertEqual(Fraction(1,2) * (-0.5), -0.25)

class TestDeterminant (unittest.TestCase):
    def testIntSize2(self):
        self.assertEqual(determinant(((1,0), (0,1))), 1)
        self.assertEqual(determinant(((1,2), (3,4))), 1*4-2*3)
        self.assertEqual(determinant(((1.0,0), (0,1))), 1.0)
        self.assertEqual(determinant(((Fraction(1,2),3), (1,5))), Fraction(-1,2))
        self.assertEqual(determinant(((Fraction(1,2),3), (1.0,5))), -0.5)

        det = determinant(((VarDbl(1,1e-3), VarDbl(2, 1e-2)), (VarDbl(3), VarDbl(4, 0.1))))
        self.assertEqual(det.value(), 1*4-2*3)
        self.assertEqual(det.variance(), (VarDbl(1,1e-3)*VarDbl(4, 0.1) - VarDbl(2, 1e-2)*VarDbl(3)).variance())

    def testIntSize3(self):
        self.assertEqual(determinant(((1,0,0), (0,1,0), (0,0,1))), 1)
        self.assertEqual(determinant(((1,2,3), (4,5,6), (7,8,9))), 
                         1*5*9 - 1*6*8 + 2*6*7 - 2*4*9 + 3*4*8 - 3*5*7)
        self.assertEqual(determinant(((1.0,0,0.0), (0,1,0.0), (0.0,0,1))), 1.0)
        self.assertEqual(determinant(((Fraction(1,2),3,4), (1,5,Fraction(-7,6)), (-8,9,6))), 
                         Fraction(1,2)*5*6 - Fraction(1,2)*Fraction(-7,6)*9 +
                         3*Fraction(-7,6)*(-8) - 3*1*6 +
                         4*1*9 - 4*5*(-8))
        self.assertEqual(determinant(((Fraction(1,2),3,4), (1,5,Fraction(-7,6)), (-8,9,6))), 
                         1/2*5*6 - 1/2*(-7/6)*9 +
                         3*(-7/6)*(-8) - 3*1*6 +
                         4*1*9 - 4*5*(-8))

        ssMatrix = ((VarDbl(1,1e-3),  VarDbl(2,1e-2), VarDbl(3,1e-4)),
                    (VarDbl(-4),     VarDbl(-5,1e-1), VarDbl(6,1e-5)),
                    (VarDbl(7,1e-7),  VarDbl(8,1e-8), VarDbl(9,1e-9)))
        det = determinant(ssMatrix)
        self.assertEqual(det.value(), 1*(-5)*9 - 1*6*8 + 2*6*7 - 2*(-4)*9 + 3*(-4)*8 - 3*(-5)*7)
        self.assertLess(det.variance(), (
                          VarDbl(1,1e-3)*VarDbl(-5,1e-1)*VarDbl(9, 1e-9) - VarDbl(1,1e-3)*VarDbl(6, 1e-5)*VarDbl(8, 1e-8) +
                          VarDbl(2, 1e-2)*VarDbl(6, 1e-5)*VarDbl(7,1e-7) - VarDbl(2,1e-2)*VarDbl(-4)*VarDbl(9,1e-9) +
                          VarDbl(3,1e-4)*VarDbl(-4)*VarDbl(8,1e-8) - VarDbl(3,1e-4)*VarDbl(-5,1e-1)*VarDbl(7,1e-7)
                         ).variance())
        self.assertAlmostEqual(det.variance(), 
                               1e-6*1e-2*1e-18 + 1e-4*1e-10*1e-14 + 1e-8*0*1e-16 + 1e-6*1e-10*1e-16 + 1e-4*0*1e-18 + 1e-8*1e-2*1e-14 +
                               1*1e-2*1e-18 + 1e-6*25*1e-18 + 1e-6*1e-2*81 +
                               4*1e-10*1e-14 + 1e-4*36*1e-14 + 1e-4*1e-10*49 +
                               9*0*1e-16 + 1e-8*16*1e-16 + 1e-8*0*64 + 
                               1*1e-10*1e-16 + 1e-6*36*1e-16 + 1e-6*1e-10*64 + 
                               4*0*1e-18 + 1e-4*16*1e-18 + 1e-4*0*81 + 
                               9*1e-2*1e-14 + 1e-8*25*1e-14 + 1e-8*1e-2*49 +
                               1e-6*(((-5)*9-6*8)**2) + 1e-2*((1*9-3*7)**2) + 1e-18*((1*(-5)-2*(-4))**2) + 
                               1e-4*(((-4)*9-6*7)**2) + 1e-10*((1*8-2*7)**2) + 1e-14*((2*6-3*(-5))**2) + 
                               1e-8*(((-4)*8-(-5)*7)**2) + 0*((2*9-3*8)**2) + 1e-16*((1*6-3*(-4))**2))


if __name__ == '__main__':
    unittest.main()