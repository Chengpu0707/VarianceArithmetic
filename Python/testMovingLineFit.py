import math
import random
import unittest

from histo import Stat
from movingLineFit import FitType, movingLineFit
from varDbl import VarDbl

class TestMovingLineFitPrecise (unittest.TestCase):

    def testFlat(self):
        sInput = [1] * 10
        sExpect = tuple([(1.0, 0.0)] * 6)
        ssFit = movingLineFit(sInput, 2, FitType.LOCAL), \
                movingLineFit(sInput, 2, FitType.MOVING_NO_VAR_ADJ), \
                movingLineFit(sInput, 2, FitType.MOVING)
        self.assertTupleEqual(ssFit[0], sExpect)
        self.assertTupleEqual(ssFit[1], ssFit[0])
        self.assertTupleEqual(ssFit[2], ssFit[0])

    def testSlope(self):
        sInput = [i for i in range(-5, 6)]
        sExpect = tuple([(float(i), 1.0) for i in range(-3,4)])
        ssFit = movingLineFit(sInput, 2, FitType.LOCAL), \
                movingLineFit(sInput, 2, FitType.MOVING_NO_VAR_ADJ), \
                movingLineFit(sInput, 2, FitType.MOVING)
        self.assertTupleEqual(ssFit[0], sExpect)
        self.assertTupleEqual(ssFit[1], ssFit[0])
        self.assertTupleEqual(ssFit[2], ssFit[0])

    def testSlopeToFlat(self):
        sInput = [-i for i in range(5)] + [-4]*5
        sExpect = tuple([(-2.0, -1.0), (-2.8, -0.8), (-3.4, -0.5), (-3.8, -0.2), (-4.0, 0.0), (-4.0, 0.0)])
        ssFit = movingLineFit(sInput, 2, FitType.LOCAL), \
                movingLineFit(sInput, 2, FitType.MOVING_NO_VAR_ADJ), \
                movingLineFit(sInput, 2, FitType.MOVING)
        self.assertEqual(len(sExpect), len(ssFit[0]))
        for i in range(len(sExpect)):
            try:
                self.assertAlmostEqual(sExpect[i][0], ssFit[0][i][0])
            except AssertionError as ex:
                raise ex
            try:
                self.assertAlmostEqual(sExpect[i][1], ssFit[0][i][1])
            except AssertionError as ex:
                raise ex

class TestMovingLineFitImprecise (unittest.TestCase):

    def verify(self, sExpect, noise, ssFit):
        self.assertEqual(len(ssFit), 3)
        for i, sFit in enumerate(ssFit):
            self.assertEqual(len(sExpect), len(sFit))
            for j, expect in enumerate(sExpect):
                try:
                    self.assertAlmostEqual(expect[0].value(), sFit[j][0].value(), delta=noise*3)
                except AssertionError as ex:
                    raise ex
                try:
                    self.assertAlmostEqual(expect[1].value(), sFit[j][1].value(), delta=noise*3)
                except AssertionError as ex:
                    raise ex
                
                match i:
                    case 0:
                        try:
                            self.assertAlmostEqual(expect[0].uncertainty(), sFit[j][0].uncertainty())
                        except AssertionError as ex:
                            raise ex
                        try:
                            self.assertAlmostEqual(expect[1].uncertainty(), sFit[j][1].uncertainty())
                        except AssertionError as ex:
                            raise ex
                    case 1:
                        try:
                            self.assertGreaterEqual(sFit[j][0].uncertainty(), ssFit[0][j][0].uncertainty())
                        except AssertionError as ex:
                            raise ex
                        try:
                            self.assertGreaterEqual(sFit[j][1].uncertainty(), ssFit[0][j][1].uncertainty())
                        except AssertionError as ex:
                            raise ex
                    case 2:
                        try:
                            self.assertAlmostEqual(sFit[j][0].uncertainty(), ssFit[0][j][0].uncertainty())
                        except AssertionError as ex:
                            raise ex
                        try:
                            self.assertAlmostEqual(sFit[j][1].uncertainty(), ssFit[0][j][1].uncertainty())
                        except AssertionError as ex:
                            raise ex
                    case _:
                        self.fail()


    def testFlat(self):
        half = 2
        noise = 1e-3
        noise0 = math.sqrt(1 /(2 * half + 1)) *noise
        noise1 = math.sqrt(3 /half /(half + 1) /(2 * half + 1)) *noise
        sInput = [VarDbl(1 + noise * random.normalvariate(), noise) for i in range(10)]
        sExpect = tuple([(VarDbl(1.0, noise0), VarDbl(0.0, noise1))] * 6)
        ssFit = movingLineFit(sInput, half, FitType.LOCAL), \
                movingLineFit(sInput, half, FitType.MOVING_NO_VAR_ADJ), \
                movingLineFit(sInput, half, FitType.MOVING)
        self.verify(sExpect, noise, ssFit)

    def testSlope(self):
        half = 2
        noise = 1e-1
        noise0 = math.sqrt(1 /(2 * half + 1)) *noise
        noise1 = math.sqrt(3 /half /(half + 1) /(2 * half + 1)) *noise
        sInput = [VarDbl(i + noise * random.normalvariate(), noise) for i in range(-5,6)]
        sExpect = tuple([(VarDbl(i, noise0), VarDbl(1, noise1)) for i in range(-3,4)])
        ssFit = movingLineFit(sInput, half, FitType.LOCAL), \
                movingLineFit(sInput, half, FitType.MOVING_NO_VAR_ADJ), \
                movingLineFit(sInput, half, FitType.MOVING)
        self.verify(sExpect, noise, ssFit)

    def testSlopeToFlat(self):
        half = 2
        noise = 1e-1
        noise0 = math.sqrt(1 /(2 * half + 1)) *noise
        noise1 = math.sqrt(3 /half /(half + 1) /(2 * half + 1)) *noise
        sInput = [VarDbl(-i + noise * random.normalvariate(), noise) for i in range(5)] + \
                 [VarDbl(-4 + noise * random.normalvariate(), noise) for i in range(5)]
        sExpect = tuple([(VarDbl(-2.0, noise0), VarDbl(-1.0, noise1)), 
                         (VarDbl(-2.8, noise0), VarDbl(-0.8, noise1)), 
                         (VarDbl(-3.4, noise0), VarDbl(-0.5, noise1)), 
                         (VarDbl(-3.8, noise0), VarDbl(-0.2, noise1)), 
                         (VarDbl(-4.0, noise0), VarDbl(0.0, noise1)), 
                         (VarDbl(-4.0, noise0), VarDbl(0.0, noise1))])
        ssFit = movingLineFit(sInput, half, FitType.LOCAL), \
                movingLineFit(sInput, half, FitType.MOVING_NO_VAR_ADJ), \
                movingLineFit(sInput, half, FitType.MOVING)
        self.verify(sExpect, noise, ssFit)

    def testDump(self):
        half = 2
        noise = 0.2
        shock = 10
        length = 10
        
        full = 2*half
        sFit0 = list(range(length - full)) + list(range(10, 5*length, -1))
        sFit1 = [1]*(length - full) + [-1]*length*4
        sStat0 = [Stat() for i in range(length*5 - full)]
        sStat1 = [Stat() for i in range(length*5 - full)]

        for repeat in range(1000):
            sInput = [VarDbl(i + noise * random.normalvariate(), noise) for i in range(length)] + \
                    [VarDbl(length - i + noise * random.normalvariate(), noise) for i in range(length)] + \
                    [VarDbl(-i + shock * noise * random.normalvariate(), noise) for i in range(length)] + \
                    [VarDbl(-length - i + noise * random.normalvariate(), noise) for i in range(length)] + \
                    [VarDbl(-length - i + noise * random.normalvariate(), noise) for i in range(length)]
            ssFit = movingLineFit(sInput, half, FitType.LOCAL), \
                    movingLineFit(sInput, half, FitType.MOVING_NO_VAR_ADJ), \
                    movingLineFit(sInput, half, FitType.MOVING)
            self.assertEqual(len(ssFit[0]), length*5 - full)
            for i, fit0 in enumerate(ssFit[0]):
                try:
                    self.assertAlmostEqual(fit0[0].value(), ssFit[2][i][0].value())
                except AssertionError as ex:
                    raise ex
                try:
                    self.assertAlmostEqual(fit0[1].value(), ssFit[2][i][1].value())
                except AssertionError as ex:
                    raise ex
            
                try:
                    self.assertAlmostEqual(fit0[0].uncertainty(), ssFit[2][i][0].uncertainty())
                except AssertionError as ex:
                    raise ex
                try:
                    self.assertAlmostEqual(fit0[1].uncertainty(), ssFit[2][i][1].uncertainty())
                except AssertionError as ex:
                    raise ex
                
                fit2 = ssFit[2][i]
                sStat0[i].accum((fit2[0].value() - sFit0[0])/fit2[0].uncertainty())
                sStat1[i].accum((fit2[1].value() - sFit1[0])/fit2[1].uncertainty())

        with open('./Python/Output/MovingLineFitStat.txt', 'w') as f:
            f.write('Time Index\tInput Value\tInput Uncertainty'
                    '\tFit 0 Deviation\tFit 0 Mean\tFit 0 Minimum\tFit 0 Maxium'
                    '\tFit 1 Deviation\tFit 1 Mean\tFit 1 Minimum\tFit 1 Maxium'
                    '\n')
            for i in range(full):
                f.write(f'{i}\t{sInput[i].value()}\t{sInput[i].uncertainty()}'
                        '\tnan\tnan\tnan\tnan'
                        '\tnan\tnan\tnan\tnan'
                        '\n')
            for i in range(len(ssFit[0])):
                try:
                    f.write(f'{i + full}\t{sInput[i + full].value()}\t{sInput[i + full].uncertainty()}')
                except BaseException as ex:
                    raise ex
                try:
                    f.write(f'\t{sStat0[i].dev()}\t{sStat0[i].mean()}\t{sStat0[i].min()}\t{sStat0[i].max()}')
                    f.write(f'\t{sStat1[i].dev()}\t{sStat1[i].mean()}\t{sStat1[i].min()}\t{sStat1[i].max()}')
                except BaseException as ex:
                    raise ex
                f.write('\n')


        with open('./Python/Output/MovingLineFit.txt', 'w') as f:
            f.write('Time Index\tInput Value\tInput Uncertainty'
                    '\tLocal 0 Value\tLocal 0 Uncertainty\tLocal 1 Value\tLocal 1 Uncertainty'
                    '\tUnadjusted 0 Value\tUnadjusted 0 Uncertainty\tUnadjusted 1 Value\tUnadjusted 1 Uncertainty'
                    '\tAdjusted 0 Value\tAdjusted 0 Uncertainty\tAdjusted 1 Value\tAdjusted 1 Uncertainty'
                    '\n')
            for i in range(full):
                f.write(f'{i}\t{sInput[i].value()}\t{sInput[i].uncertainty()}'
                        '\tnan\tnan\tnan\tnan'
                        '\tnan\tnan\tnan\tnan'
                        '\tnan\tnan\tnan\tnan'
                        '\n')
            for i in range(len(ssFit[0])):
                try:
                    f.write(f'{i + full}\t{sInput[i + full].value()}\t{sInput[i + full].uncertainty()}')
                except BaseException as ex:
                    raise ex
                for j in range(3):
                    for k in range(2):
                        try:
                            f.write(f'\t{ssFit[j][i][k].value()}\t{ssFit[j][i][k].uncertainty()}')
                        except BaseException as ex:
                            raise ex
                f.write('\n')




        


if __name__ == '__main__':
    unittest.main()


