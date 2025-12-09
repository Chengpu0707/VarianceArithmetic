
import math
import logging
import os
import unittest

from taylor import Taylor, NotMonotonicException, NotFiniteException
from varDbl import VarDbl

logger = logging.getLogger(__name__)


class TestPolynominial (unittest.TestCase):

    def test_poly_0(self):
        res = Taylor.polynominal1d(VarDbl(), (1,))
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

        res = Taylor.polynominal1d(VarDbl(1), (0,))
        self.assertEqual(res.value(), 0)
        self.assertEqual(res.uncertainty(), 0)

        res = Taylor.polynominal1d(VarDbl(0, 0.1), (1,))
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

        res = Taylor.polynominal1d(VarDbl(0, 0.1), (VarDbl(0, 0.1),))
        self.assertEqual(res.value(), 0)
        self.assertAlmostEqual(res.uncertainty(), 0.1)

        res = Taylor.polynominal1d(VarDbl(1), (1,))
        self.assertEqual(res.value(), 1)
        self.assertEqual(res.uncertainty(), 0)

        res = Taylor.polynominal1d(VarDbl(1), (-2,))
        self.assertEqual(res.value(), -2)
        self.assertEqual(res.uncertainty(), 0)

        res = Taylor.polynominal1d(VarDbl(1, 0.1), (-2,))
        self.assertEqual(res.value(), -2)
        self.assertAlmostEqual(res.uncertainty(), 0)

        res = Taylor.polynominal1d(VarDbl(2), (VarDbl(1, 0.1),))
        self.assertEqual(res.value(), 1)
        self.assertAlmostEqual(res.uncertainty(), 0.1)

    def test_poly_1(self):
        res = Taylor.polynominal1d(VarDbl(0, 1/8), (0,1))
        self.assertEqual(res.value(), 0)
        self.assertAlmostEqual(res.uncertainty(), 1/8, delta=1e-6)

        res = Taylor.polynominal1d(VarDbl(0, 1/8), (1,1))
        self.assertEqual(res.value(), 1)
        self.assertAlmostEqual(res.uncertainty(), 1/8, delta=1e-6)

        res = Taylor.polynominal1d(VarDbl(-2, 1/8), (1,1))
        self.assertEqual(res.value(), -1)
        self.assertAlmostEqual(res.uncertainty(), 1/8, delta=1e-6)

    def test_poly_2(self):
        for value in (0, 1, -1, 2, -2, 0.25, -0.25):
            try:
                res = Taylor.polynominal1d(VarDbl(value, 0.5), (0,0,1))
                res2 = VarDbl(value**2 + 0.5**2, math.sqrt(4 *value**2 *0.5**2 + 2 *0.5**4))
                self.assertAlmostEqual(res.value(), res2.value(), delta=1e-5)
                self.assertAlmostEqual(res.uncertainty(), res2.uncertainty(), delta=5e-5)
            except BaseException as ex:
                raise ex
            try:
                res = Taylor.polynominal1d(VarDbl(value, 0.5), (1,2,1))
                res2 = VarDbl((value + 1)**2 + 0.5**2, math.sqrt(4 *(value + 1)**2 *0.5**2 + 2 *0.5**4))
                self.assertAlmostEqual(res.value(), res2.value(), delta=1e-5)
                self.assertAlmostEqual(res.uncertainty(), res2.uncertainty(), delta=5e-5)
            except BaseException as ex:
                raise ex
            try:
                res = Taylor.polynominal1d(VarDbl(value, 0.5), (1,-2,1))
                res2 = VarDbl((value - 1)**2 + 0.5**2, math.sqrt(4 *(value - 1)**2 *0.5**2 + 2 *0.5**4))
                self.assertAlmostEqual(res.value(), res2.value(), delta=1e-5)
                self.assertAlmostEqual(res.uncertainty(), res2.uncertainty(), delta=5e-5)
            except BaseException as ex:
                raise ex

    def test_poly_3(self):
        for value in (0, 1, -1, 2, -2, 0.25, -0.25):
            try:
                res = Taylor.polynominal1d(VarDbl(value, 0.5), (0,0,0,1))
                res2 = VarDbl(value, 0.5)**3
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            try:
                res = Taylor.polynominal1d(VarDbl(value, 0.5), (1,3,3,1))
                res2 = VarDbl(1 + value, 0.5)**3
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            try:
                res = Taylor.polynominal1d(VarDbl(value, 0.5), (1,-3,3,-1))
                res2 = VarDbl(1 - value, 0.5)**3
                self.assertEqual(res.value(), res2.value())
                self.assertEqual(res.uncertainty(), res2.uncertainty())
            except BaseException as ex:
                raise ex
            


class TestPolyOne (unittest.TestCase):
    MAX_ORDER = 224
    sPoly = [1] * MAX_ORDER

    def test_precise_x(self):
        inv = Taylor.pow(1 - VarDbl(0.7), -1, dumpPath='./Python/Output/Pow_0.3_-1.txt')
        self.assertAlmostEqual(inv.value(), 10/3)
        self.assertAlmostEqual(inv.uncertainty() / 7.569487334792672e-16, 1)
        res = Taylor.polynominal1d(VarDbl(0.7), TestPolyOne.sPoly, dumpPath='./Python/Output/PolyOne_0.7.txt')
        self.assertAlmostEqual(res.value(), 10/3)
        self.assertAlmostEqual(res.uncertainty(), 7.569487334792672e-16)
        self.assertAlmostEqual(res.value(), inv.value(), delta=res.uncertainty()*10)

        inv = Taylor.pow(1 - VarDbl(0.8), -1, dumpPath='./Python/Output/Pow_0.2_-1.txt')
        self.assertAlmostEqual(inv.value(), 5)
        self.assertAlmostEqual(inv.uncertainty() / 1.682504257894355e-15, 1)
        res = Taylor.polynominal1d(VarDbl(0.8), TestPolyOne.sPoly, dumpPath='./Python/Output/PolyOne_0.8.txt')
        self.assertAlmostEqual(res.value(), 5)
        self.assertAlmostEqual(res.uncertainty(), 1.6825042578220935e-15)
        self.assertAlmostEqual(res.value(), inv.value(), delta=res.uncertainty())

        inv = Taylor.pow(1 - VarDbl(0.99), -1, dumpPath='./Python/Output/Pow_0.01_-1.txt')
        self.assertAlmostEqual(inv.value(), 99.99999999999991)
        self.assertAlmostEqual(inv.uncertainty(), 6.410351214808837e-13)
        res = Taylor.polynominal1d(VarDbl(0.99), TestPolyOne.sPoly, dumpPath='./Python/Output/PolyOne_0.99.txt')
        self.assertAlmostEqual(res.value(), 89.47350981516396)
        self.assertAlmostEqual(res.uncertainty(), 4.209234934156554e-13)

        with open('./Python/Output/PolyOne_x.txt', 'w') as f:
            f.write('x\tInverse Value\tInverse Uncertainty\tPoly Value\tPoly Uncertainty\tValue Diff\tValue LSV\tLSV Diff\tOrder\n')
            for i in range(-99, 99):
                x = i / 100
                inv = Taylor.pow(1 - VarDbl(x), -1)
                poly = Taylor.polynominal1d(VarDbl(x), TestPolyOne.sPoly)
                diff = poly.value() - inv.value()
                lsv = math.ulp(inv.value())
                exp = 1
                j = 0
                while exp >= lsv:
                    exp *= abs(x)
                    j += 1
                f.write(f'{x}\t{inv.value()}\t{inv.uncertainty()}\t{poly.value()}\t{poly.uncertainty()}'
                        f'\t{diff}\t{lsv}\t{diff / lsv}\t{j}\n')

    def test_precise_order(self):
        inv = Taylor.pow(1 - VarDbl(0.6), -1)
        self.assertAlmostEqual(inv.value(), 2.5)
        self.assertAlmostEqual(inv.uncertainty() / 4.00614133453171e-16, 1)

        res = Taylor.polynominal1d(VarDbl(0.6), [1]*69)
        self.assertAlmostEqual((res.value() - inv.value()) / math.ulp(inv.value()), -1)
        self.assertAlmostEqual(res.uncertainty() / 4.756365132018834e-16, 1)

        res = Taylor.polynominal1d(VarDbl(0.6), [1]*70)
        self.assertAlmostEqual((res.value() - inv.value()) / math.ulp(inv.value()), 0)
        self.assertAlmostEqual(res.uncertainty() / 4.0061413345316546e-16, 1)

        res = Taylor.polynominal1d(VarDbl(0.6), [1]*71)
        self.assertAlmostEqual((res.value() - inv.value()) / math.ulp(inv.value()), 1)
        self.assertAlmostEqual(res.uncertainty() / 4.756365132018883e-16, 1)

        with open('./Python/Output/PolyOne_order.txt', 'w') as fw:
            fw.write('x\tOrder\tPoly Value\tPoly Uncertainty\tValue Diff\tPoly LSV\n')
            for x in (-0.7, -0.6, 0.6, 0.7):
                inv = Taylor.pow(1 - VarDbl(x), -1)
                for n in range(2, TestPolyOne.MAX_ORDER):
                    poly = Taylor.polynominal1d(VarDbl(x), [1]*n)
                    fw.write(f'{x}\t{n}\t{poly.value()}\t{poly.uncertainty()}\t{poly.value() - inv.value()}\t{math.ulp(poly.value())}\n')

    def test_imprecise(self):
        with self.assertRaises(NotFiniteException):
            Taylor.polynominal1d(VarDbl(0.7, 0.37), TestPolyOne.sPoly)
        res = Taylor.polynominal1d(VarDbl(0.7, 0.36), TestPolyOne.sPoly)
        self.assertAlmostEqual(res.value(), 5.00e+81, delta=1e79)
        self.assertAlmostEqual(res.uncertainty(), 1.45e+85, delta=1e83)

        with self.assertRaises(NotMonotonicException):
            Taylor.pow(1 - VarDbl(0.7, 0.061), -1)
        res = Taylor.pow(1 - VarDbl(0.7, 0.060), -1)
        self.assertAlmostEqual(res.value(), 10/3 + 0.1541670)
        self.assertAlmostEqual(res.uncertainty(), 0.8494772)
        res = Taylor.polynominal1d(VarDbl(0.7, 0.060), TestPolyOne.sPoly)
        self.assertAlmostEqual(res.value(), 10/3 + 0.1541102)
        self.assertAlmostEqual(res.uncertainty(), 0.8299787)

        with open('./Python/Output/PolyOne_x_0.05.txt', 'w') as f:
            f.write('x\tInverse Value\tInverse Uncertainty\tPoly Value\tPoly Uncertainty\tValue Diff\n')
            for i in range(-75, 75):
                x = i / 100
                inv = Taylor.pow(1 - VarDbl(x, 0.05), -1)
                poly = Taylor.polynominal1d(VarDbl(x, 0.05), TestPolyOne.sPoly)
                f.write(f'{x}\t{inv.value()}\t{inv.uncertainty()}\t{poly.value()}\t{poly.uncertainty()}\t{poly.value() - inv.value()}\n')



if __name__ == '__main__':
    unittest.main()

