from fractions import Fraction
import itertools
import numpy
import unittest

from matrix import permutSign, isSquareMatrix, createIntMatrix, addNoise, linear, multiply, adjugate
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
        Can not use numpy to hold int matrix and to retain int for determinant
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

class TestMultiply (unittest.TestCase):
    def testIntSize2(self):
        ssId = ((1,0),(0,1))
        ssTr = ((0,1),(1,0))
        ssMat = ((1,-2),(-3,4))

        self.assertTupleEqual(multiply(ssId, ssId), ssId)

        self.assertTupleEqual(multiply(ssTr, ssId), ssTr)
        self.assertTupleEqual(multiply(ssId, ssTr), ssTr)
        self.assertTupleEqual(multiply(ssTr, ssTr), ssId)

        self.assertTupleEqual(multiply(ssId, ssMat), ssMat)
        self.assertTupleEqual(multiply(ssMat, ssId), ssMat)

        self.assertTupleEqual(multiply(ssTr, ssMat), ((-3,4),(1,-2)))
        self.assertTupleEqual(multiply(ssMat, ssTr), ((-2,1),(4,-3)))

        self.assertTupleEqual(multiply(ssMat, ssMat), 
                              ((1+6, -2-8), (-3-12, 6+16)))
        
    def testDiffSize(self):
        with self.assertRaises(ValueError):
            multiply(((1,0),(0,1)), ((1,0,0),(0,10),(0,0,1)))

class TestLinear (unittest.TestCase):
    def verifyValue(self, val, ret):
        try:
            if eType := type(val) == VarDbl:
                self.assertAlmostEqual(val.value(), ret.value())
                self.assertAlmostEqual(val.variance(), ret.variance())
            elif eType == float:
                self.assertAlmostEqual(val, ret)
            else:
                self.assertEqual(val, ret)
        except AssertionError as ex:
            raise ex

    def testIntSize2(self):
        ssMat = ((1,-2),(-3,4))
        self.assertTupleEqual(linear(ssMat), ssMat)
        self.assertTupleEqual(linear(ssMat, 2), ((2,-4),(-6,8)))
        self.assertTupleEqual(linear(ssMat, 2, -1), ((1,-5),(-7,7)))

    def testVardblSize2(self):
        ssMat = (VarDbl(1,1e-3), VarDbl(2, 1e-2)), (VarDbl(3), VarDbl(4, 0.1))
        ssRes = linear(ssMat, -2, 1)
        self.verifyValue(ssRes[0][0], VarDbl(-1, 2e-3))
        self.verifyValue(ssRes[0][1], VarDbl(-3, 2e-2))
        self.verifyValue(ssRes[1][0], VarDbl(-5))
        self.verifyValue(ssRes[1][1], VarDbl(-7, 0.2))
 

class TestAdjugate (unittest.TestCase):
        
    def verifyValue(self, val, ret):
        try:
            if eType := type(val) == VarDbl:
                self.assertAlmostEqual(val.value(), ret.value())
                self.assertAlmostEqual(val.variance(), ret.variance())
            elif eType == float:
                self.assertAlmostEqual(val, ret)
            else:
                self.assertEqual(val, ret)
        except AssertionError as ex:
            raise ex

    def verifyIdentity(self, det, ssId):
        size = len(ssId)
        for i in range(size):
            for j in range(size):
                if eType := type(ssId[i][j]) == VarDbl:
                    self.assertAlmostEqual(ssId[i][j].value(), det.value() if i == j else 0)
                    #self.assertAlmostEqual(ssId[i][j].variance(), det.variance() if i == j else 0)
                elif eType == float:
                    self.assertAlmostEqual(ssId[i][j], det if i == j else 0)
                else:
                    self.assertEqual(ssId[i][j], det if i == j else 0)

    def verify(self, ssMat, det, ssAdj):
        size = len(ssMat)
        self.verifyIdentity(det, multiply(ssMat, ssAdj))
        ret, ssRet = adjugate(ssMat)
        self.verifyValue(det, ret)
        self.verifyIdentity(det, multiply(ssMat, ssRet))
        for i in range(size):
            for j in range(size):
                self.verifyValue(ssAdj[i][j], ssRet[i][j])


    def testIntSize2(self):
        self.verify(((1,0), (0,1)), 1, ((1,0), (0,1)))
        self.verify(((0,1), (1,0)), -1, ((0,-1), (-1,0)))
        self.verify(((1,2), (3,4)), -2, ((4,-2), (-3,1)))

    def testFloatSize2(self):
        self.verify(((1,2.0), (3,4)), -2.0, ((4,-2.0), (-3,1)))

    def testIntFractionSize2(self):
        self.verify(((1,2), (3,Fraction(1/2))), Fraction(-11, 2), ((Fraction(1/2),-2), (-3,1)))

    def testIntFractionFloatSize2(self):
        self.verify(((1,2.0), (3,Fraction(1/2))), -11/2, ((Fraction(1/2),-2.0), (-3,1)))

    def testVarDblSize2(self):
        ssMat = ((VarDbl(1,1e-3), VarDbl(2, 1e-2)), 
                 (VarDbl(3), VarDbl(4, 0.1)))
        det = VarDbl(1,1e-3)*VarDbl(4, 0.1) - VarDbl(2, 1e-2)*VarDbl(3)
        ssAdj = ((VarDbl(4, 4e-3), -2), 
                 (VarDbl(-3, 3e-2), VarDbl(1,0.1)))
        self.verify(ssMat, det, ssAdj)

    def testIntSize3(self):
        self.verify(((1,0,0), (0,1,0), (0,0,1)), 1, ((1,0,0), (0,1,0), (0,0,1)))
        self.verify(((1,2,3), (-4,-5,6), (7,8,9)), 72, ((-93,6,27), (78,-12,-18), (3,6,3)))

    def testVarDblSize3_noTracing(self):
        noTracing = VarDbl(1,1e-3)*VarDbl(-5,1e-1)*VarDbl(9,1e-9) - VarDbl(1,1e-3)*VarDbl(6,1e-5)*VarDbl(8,1e-8) +\
                    VarDbl(2,1e-2)*VarDbl(6,1e-5)*VarDbl(7,1e-7) - VarDbl(2,1e-2)*VarDbl(-4)*VarDbl(9,1e-9) +\
                    VarDbl(3,1e-4)*VarDbl(-4)*VarDbl(8,1e-8) - VarDbl(3,1e-4)*VarDbl(-5,1e-1)*VarDbl(7,1e-7)
        self.verifyValue(noTracing, VarDbl(72, 5.530352330904206, True))

    def testVarDblSize3(self):
        ssMat3 = ((VarDbl( 1,1e-3), VarDbl( 2,1e-2), VarDbl(3,1e-4)),
                  (-4,              VarDbl(-5,1e-1), VarDbl(6,1e-5)),
                  (VarDbl( 7,1e-7), VarDbl( 8,1e-8), VarDbl(9,1e-9)))
        val3 = 1*(-5)*9 - 1*6*8 + 2*6*7 - 2*(-4)*9 + 3*(-4)*8 - 3*(-5)*7
        self.assertEqual(val3, 72)
        var3 = 1e-6*1e-2*1e-18 + 1e-4*1e-10*1e-14 + 1e-8*0*1e-16 + 1e-6*1e-10*1e-16 + 1e-4*0*1e-18 + 1e-8*1e-2*1e-14 +\
                1*1e-2*1e-18 + 1e-6*25*1e-18 + 1e-6*1e-2*81 +\
                4*1e-10*1e-14 + 1e-4*36*1e-14 + 1e-4*1e-10*49 +\
                9*0*1e-16 + 1e-8*16*1e-16 + 1e-8*0*64 +\
                1*1e-10*1e-16 + 1e-6*36*1e-16 + 1e-6*1e-10*64 +\
                4*0*1e-18 + 1e-4*16*1e-18 + 1e-4*0*81 +\
                9*1e-2*1e-14 + 1e-8*25*1e-14 + 1e-8*1e-2*49 +\
                1e-6*(((-5)*9-6*8)**2) + 1e-2*((1*9-3*7)**2) + 1e-18*((1*(-5)-2*(-4))**2) +\
                1e-4*(((-4)*9-6*7)**2) + 1e-10*((1*8-2*7)**2) + 1e-14*((2*6-3*(-5))**2) +\
                1e-8*(((-4)*8-(-5)*7)**2) + 0*((2*9-3*8)**2) + 1e-16*((1*6-3*(-4))**2)
        self.assertEqual(var3, 2.0570499085078198)

        ssAdj = (
            (
                VarDbl(-93, (93e-3)**2 + (1e-3*1e-1*1e-9)**2 + (1e-3*5*1e-9)**2 + (1e-3*1e-1*9)**2
                                       + (1e-3*1e-5*1e-8)**2 + (1e-3*6*1e-8)**2 + (1e-3*1e-5*8)**2, True),
                6, 
                VarDbl(27,  (27e-7)**2 + (1e-7*1e-2*1e-5)**2 + (1e-7*2*1e-5)**2 + (1e-7*1e-2*6)**2
                                       + (1e-7*1e-1*1e-4)**2 + (1e-7*5*1e-4)**2 + (1e-7*1e-1*3)**2, True),
            ), 
            (
                VarDbl(78,  (78e-2)**2 + (1e-2*1e-5*1e-7)**2 + (1e-2*6*1e-7)**2 + (1e-2*1e-5*7)**2
                                       + (1e-2*0*1e-9)**2    + (1e-2*4*1e-9)**2 + (1e-2*0*9)**2, True ),
                VarDbl(-12,  (12e-1)**2 + (1e-1*1e-3*1e-9)**2 + (1e-1*1*1e-9)**2 + (1e-1*1e-3*9)**2
                                       + (1e-1*1e-4*1e-7)**2 + (1e-1*3*1e-7)**2 + (1e-1*1e-4*7)**2, True),
                VarDbl(-18, (18e-8)**2 + (1e-8*1e-3*1e-5)**2 + (1e-8*1*1e-5)**2 + (1e-8*1e-3*6)**2, True),
            ),
            (
                VarDbl(3,   (3e-4)**2  + (1e-4*0*1e-8)**2    + (1e-4*4*1e-8)**2 + (1e-4*0*8)**2 
                                       + (1e-4*1e-1*1e-7)**2 + (1e-4*5*1e-7)**2 + (1e-4*1e-1*7)**2, True),
                VarDbl(6,   (6e-5)**2  + (1e-5*1e-2*1e-7)**2 + (1e-5*2*1e-7)**2 + (1e-5*1e-2*7)**2
                                       + (1e-5*1e-3*1e-8)**2 + (1e-5*1*1e-8)**2 + (1e-5*1e-3*8)**2, True ),
                VarDbl(3,   (3e-9)**2  + (1e-9*1e-3*1e-1)**2 + (1e-9*1*1e-1)**2 + (1e-9*1e-3*5)**2
                                       + (1e-9*1e-2*0)**2    + (1e-9*2*0)**2    + (1e-9*1e-2*4)**2, True),
            ))
        self.verify(ssMat3, VarDbl(val3, var3, True), ssAdj)


if __name__ == '__main__':
    unittest.main()