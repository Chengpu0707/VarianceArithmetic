import datetime
from fractions import Fraction
import itertools
import math
import numpy
import unittest

from histo import Histo, Stat
from matrix import ElementType, permutSign, isSquareMatrix, createIntMatrix, addNoise, linear, multiply, adjugate
from taylor import LossUncertaintyException
from varDbl import VarDbl, UncertaintyException


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

    def testAddNoise(self):
        for size in range(4, Adjugate.MAX_SIZE):
            noise = 1e-3
            ssOrg = createIntMatrix(size)
            ssMat = addNoise(ssOrg, 1e-3)
            self.assertTrue(isSquareMatrix(ssMat, sType=(VarDbl,)))
            uncStat = Stat()
            valStat = Stat()
            errStat = Stat()
            nrmStat = Stat()
            try:
                for i in range(size):
                    for j in range(size):
                        valStat.accum(ssOrg[i][j])
                        self.assertEqual(type(ssMat[i][j]), VarDbl)
                        self.assertEqual(type(ssMat[i][j].value()), float)
                        unc = ssMat[i][j].uncertainty()
                        uncStat.accum(unc)
                        errStat.accum(ssMat[i][j].value() - ssOrg[i][j])
                        nrmStat.accum((ssMat[i][j].value() - ssOrg[i][j])/unc)
                dev =  ELEMENT_RANGE/math.sqrt(3)
                self.assertTrue(abs(valStat.mean()) < dev/2, f'size {size}: val_mean = {valStat.mean()} > {dev/2}, val_dev = {valStat.dev()}')
                self.assertAlmostEqual(valStat.dev(), dev, delta=dev/4, msg=f'size {size}: val_mean = {valStat.mean()}, val_dev = {valStat.dev()} vs {dev}')
                dev *= noise
                self.assertEqual(uncStat.dev(), 0)
                self.assertAlmostEqual(uncStat.mean(), dev, msg=f'size {size}: uncertainty_mean = {uncStat.mean()}, uncertainty_dev = {uncStat.dev()} vs {dev}')
                self.assertAlmostEqual(errStat.dev(), dev, delta=dev/2, msg=f'size {size}: error_mean = {errStat.mean()} > {dev/2}, error_dev = {errStat.dev()}')
                self.assertTrue(0.75 <= nrmStat.dev() <= 1.25, msg=f'size {size}: normalized_mean = {nrmStat.mean()} not in [0.75, 1.25], normalized_dev = {nrmStat.dev()}')
            except BaseException as ex:
                raise ex


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


NOISES = tuple([math.pow(10, nexp) for nexp in range(-17, 0)])
HIST_DIVIDES = 5
HIST_RANGE = 3


class Adjugate:
    '''
    A class to hold precise matrix for adjugate matrix.
    A precise matrix has all int elements.
    ELEMENT_RANGE and ROUND_RANGE so that the calculated variance does not overflow
    '''
    MAX_SIZE = 11
    ELEMENT_RANGE = 1 << 8

    @staticmethod
    def noise(noiseLevel:float) ->float:
        return Adjugate.ELEMENT_RANGE/math.sqrt(3) * noiseLevel

    sLossVar = {(size, noise): 0 for size in range(MAX_SIZE) for noise in (list(NOISES) + [0])}
    sAdjStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in (list(NOISES) + [0])}
    sAdjHisto = {(size, noise): Histo(HIST_DIVIDES, HIST_RANGE) 
                 for size in range(MAX_SIZE) for noise in (list(NOISES) + [0])}
    sFwdStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in (list(NOISES) + [0])}
    sRndStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in (list(NOISES) + [0])}

    __slots__ = ('size', 'randRange', 'ssOrg', 'ssAdj', 'detAdj')

    def __init__(self, size:int, randRange=ELEMENT_RANGE) -> None:
        self.size = size
        self.randRange = randRange
        self.ssOrg = createIntMatrix(size, randRange=randRange)
        self.detAdj, self.ssAdj = adjugate(self.ssOrg)
    
        ssId = multiply(self.ssOrg, self.ssAdj)
        detRnd, ssRnd = adjugate(self.ssAdj)
        ssId2 = multiply(self.ssAdj, ssRnd)
        for i in range(size):
            for j in range(size):
                try:
                    assert ssId[i][j] == (self.detAdj if i == j else 0)
                    assert ssId2[i][j] == (detRnd if i == j else 0)
                    assert (ssRnd[i][j] * self.detAdj) == (self.ssOrg[i][j] * detRnd)
                except AssertionError as ex:
                    raise ex



class TestAdjugate (unittest.TestCase):
        
    def verifyValue(self, val, ret):
        try:
            eType = type(val)
            if eType == VarDbl:
                self.assertAlmostEqual(val.value(), ret.value())
                self.assertAlmostEqual(val.variance(), ret.variance(), delta=val.variance()*1e-2)
            elif eType == float:
                self.assertAlmostEqual(val, ret)
            else:
                self.assertEqual(val, ret)
        except AssertionError as ex:
            raise ex

    def verifyIdentity(self, det, ssId, places=5):
        size = len(ssId)
        for i in range(size):
            for j in range(size):
                try:
                    eType = type(ssId[i][j])
                    if eType == VarDbl:
                        self.assertAlmostEqual(ssId[i][j].value(), det.value() if i == j else 0,
                                               delta=abs(det.value())*1e-6)
                    elif eType == float:
                        self.assertAlmostEqual(ssId[i][j], det if i == j else 0, places=places)
                    else:
                        self.assertEqual(ssId[i][j], det if i == j else 0)
                except AssertionError as ex:
                    raise ex

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
                 (3,              VarDbl(4, 1e-1)))
        det = VarDbl(1,1e-3)*VarDbl(4, 0.1) - VarDbl(2, 1e-2)*VarDbl(3)
        ssAdj = ((VarDbl(4,1e-1),  VarDbl(-2,1e-2)), 
                 (-3,              VarDbl( 1,1e-3)))
        self.verify(ssMat, det, ssAdj)


    def testIntSize3(self):
        self.verify(((1,0,0), (0,1,0), (0,0,1)), 1, ((1,0,0), (0,1,0), (0,0,1)))
        self.verify(((1,2,3), (-4,-5,6), (7,8,9)), 72, ((-93,6,27), (78,-12,-18), (3,6,3)))

    def testVarDblSize3_noTracing(self):
        noTracing = VarDbl( 1,1e-1)*VarDbl(-5,1e-5)*VarDbl( 9,1e-9) - VarDbl( 1,1e-1)*VarDbl( 6,1e-6)*VarDbl( 8,1e-8) +\
                    VarDbl( 2,1e-2)*VarDbl( 6,1e-6)*VarDbl( 7,1e-7) - VarDbl( 2,1e-2)*VarDbl(-4,1e-4)*VarDbl( 9,1e-9) +\
                    VarDbl( 3,1e-3)*VarDbl(-4,1e-4)*VarDbl( 8,1e-8) - VarDbl( 3,1e-3)*VarDbl(-5,1e-5)*VarDbl( 7,1e-7)
        self.verifyValue(noTracing, VarDbl(72, 43.59825805262702, True))

    def testVarDblSize3(self):
        ssMat3 = ((VarDbl( 1,1e-1), VarDbl( 2,1e-2), VarDbl( 3,1e-3)),
                  (VarDbl(-4,1e-4), VarDbl(-5,1e-5), VarDbl( 6,1e-6)),
                  (VarDbl( 7,1e-7), VarDbl( 8,1e-8), VarDbl( 9,1e-9)))
        val3 = 72
        var3 =  (1e-1*1e-5*1e-9)**2 + (1e-1*1e-6*1e-8)**2 + (1e-2*1e-6*1e-7)**2 + (1e-2*1e-4*1e-9)**2 + (1e-3*1e-4*1e-8)**2 + (1e-3*1e-5*1e-7) +\
                (1**2)*1e-5*1e-9 + 1e-1*(5**2)*1e-9 + 1e-1*1e-5*(9**2) + (1**2)*1e-6*1e-8 + 1e-1*(6**2)*1e-8 + 1e-1*1e-6*(8**2) +\
                (2**2)*1e-6*1e-7 + 1e-2*(6**2)*1e-7 + 1e-2*1e-6*(7**2) + (2**2)*1e-4*1e-9 + 1e-2*(4**2)*1e-9 + 1e-2*1e-4*(9**2) +\
                (3**2)*1e-4*1e-8 + 1e-3*(4**2)*1e-8 + 1e-3*1e-4*(8**2) + (3**2)*1e-5*1e-7 + 1e-3*(5**2)*1e-7 + 1e-3*1e-5*(7**2) +\
                (1e-1*((-5)*9-6*8))**2 + (1e-2*((-4)*9-6*7))**2 + (1e-3*((-4)*8-(-5)*7))**2 +\
                (1e-4*(2*9-3*8))**2    + (1e-5*(1*9-3*7))**2    + (1e-6*(1*8-2*7))**2 +\
                (1e-7*(2*6-3*(-5)))**2 + (1e-8*(1*6-3*(-4)))**2 + (1e-9*(1*(-5)-2*(-4)))**2
        self.assertEqual(var3, 87.09858523178215)

        ssAdj = (
            (
                VarDbl(-93, (1e-5*1e-9)**2 + (5*1e-9)**2 + (1e-5*9)**2
                          + (1e-6*1e-8)**2 + (6*1e-8)**2 + (1e-6*8)**2, True),
                VarDbl(6,   (1e-2*1e-9)**2 + (2*1e-9)**2 + (1e-2*9)**2
                          + (1e-8*1e-3)**2 + (8*1e-3)**2 + (1e-8*3)**2, True),
                VarDbl(27,  (1e-2*1e-6)**2 + (2*1e-5)**2 + (1e-2*6)**2
                          + (1e-5*1e-3)**2 + (5*1e-3)**2 + (1e-5*3)**2, True),
            ), 
            (
                VarDbl(78,  (1e-5*1e-7)**2 + (6*1e-7)**2 + (1e-5*7)**2
                          + (1e-4*1e-9)**2 + (4*1e-9)**2 + (1e-4*9)**2, True ),
                VarDbl(-12, (1e-1*1e-9)**2 + (1*1e-9)**2 + (1e-1*9)**2
                          + (1e-3*1e-7)**2 + (3*1e-7)**2 + (1e-3*7)**2, True),
                VarDbl(-18, (1e-1*1e-6)**2 + (1*1e-6)**2 + (1e-1*6)**2
                          + (1e-4*1e-3)**2 + (4*1e-3)**2 + (1e-4*3)**2, True),
            ),
            (
                VarDbl(3,   (1e-4*1e-8)**2 + (4*1e-8)**2 + (1e-4*8)**2 
                          + (1e-5*1e-7)**2 + (5*1e-7)**2 + (1e-5*7)**2, True),
                VarDbl(6,   (1e-2*1e-7)**2 + (2*1e-7)**2 + (1e-2*7)**2
                          + (1e-1*1e-8)**2 + (1*1e-8)**2 + (1e-1*8)**2, True ),
                VarDbl(3,   (1e-1*1e-5)**2 + (1*1e-5)**2 + (1e-1*5)**2
                          + (1e-2*1e-4)**2 + (2*1e-4)**2 + (1e-2*4)**2, True),
            ))
        self.verify(ssMat3, VarDbl(val3, var3, True), ssAdj)


    def roundtrip(self, adj:Adjugate, noise:float, ssOrg:tuple[tuple[ElementType]],
                  verify:bool=False):
        detAdj, ssAdj = adjugate(ssOrg)
        ssId = multiply(ssOrg, ssAdj)
        detRnd, ssRnd = adjugate(ssAdj)
            
        if verify:
            self.verifyIdentity(detAdj, multiply(ssOrg, ssAdj))
            self.verifyIdentity(detRnd, multiply(ssAdj, ssRnd), places=5)

        def accum(value, expected, adjStat:Stat, adjHisto:Histo):
            if type(value) == VarDbl:
                adjUnc = value.uncertainty()
                adjStat.accum(adjUnc)
                adjHisto.accum((value.value() - expected)/adjUnc)
                return True
            else:
                Adjugate.sLossVar[(adj.size, noise)] += 1
                return False

        adjStat = Adjugate.sAdjStat[(adj.size, noise)]
        adjHisto = Adjugate.sAdjHisto[(adj.size, noise)]
        accum(detAdj, adj.detAdj, adjStat, adjHisto)

        precStat = Stat()
        adjStat = Adjugate.sAdjStat[(adj.size, noise)]
        adjHisto = Adjugate.sAdjHisto[(adj.size, noise)]
        fwdStat = Adjugate.sFwdStat[(adj.size, noise)]
        rndStat = Adjugate.sRndStat[(adj.size, noise)]
        for i in range(adj.size):
            for j in range(adj.size):
                if accum(ssAdj[i][j], adj.ssAdj[i][j], adjStat, adjHisto) and ssAdj[i][j].value():
                    precStat.accum(ssAdj[i][j].uncertainty() / abs(ssAdj[i][j].value()))

                diff = ssId[i][j] - detAdj if i == j else ssId[i][j]
                if type(diff) == VarDbl:
                    fwdStat.accum(diff.value()/ diff.uncertainty())
                else:
                    Adjugate.sLossVar[(adj.size, noise)] += 1

                diff = detAdj * ssRnd[i][j] - detRnd * ssOrg[i][j]
                if type(diff) == VarDbl:
                    rndStat.accum(diff.value()/diff.uncertainty())
                else:
                    Adjugate.sLossVar[(adj.size, noise)] += 1

        return detAdj, detRnd, precStat


    def testRoundtrip(self):
        for size in range(4, 7):
            now = datetime.datetime.now()
            adj = Adjugate(size)
            for noise in (1e-3, 1e-9):
                for repeat in range(Adjugate.MAX_SIZE - size):
                    ssVarOrg = addNoise(adj.ssOrg, Adjugate.noise(noise))
                    self.roundtrip(adj, noise, ssVarOrg, verify=True)

                adjHisto = Adjugate.sAdjHisto[(adj.size, noise)]
                rndStat = Adjugate.sRndStat[(adj.size, noise)]
                try:
                    self.assertGreater(adjHisto.stat().dev(), 0.5)
                    self.assertLess(adjHisto.stat().dev(), 2)
                    self.assertLess(abs(adjHisto.stat().mean()), 0.3)
                    self.assertLess(rndStat.dev(), 0.1)
                    self.assertLess(rndStat.mean(), 0.1)
                except BaseException as ex:
                    raise ex
                
            adj = Adjugate(size, randRange=(1 << 14))
            for repeat in range(Adjugate.MAX_SIZE - size):
                ssVarOrg = addNoise(adj.ssOrg, 0)
                self.roundtrip(adj, 0, ssVarOrg)
            adjHisto = Adjugate.sAdjHisto[(adj.size, 0)]
            if adjHisto.stat().dev():
                try:
                    self.assertGreater(adjHisto.stat().dev(), 0.1)
                    self.assertLess(adjHisto.stat().mean(), 10)
                except BaseException as ex:
                    raise ex
            print(f'{size}: {(datetime.datetime.now() - now).total_seconds()}')

    def testException(self):
        '''
        The adjugate matrix has order of magnitude of 2^{18*6} == 10^{31}.  
        The actual measured is 10^{28}~10^{26}.
        The roundtrip matrix has order of magnitude of 10^{28*6} == 10^{206}.  
        The square of the value is 10^{412}. The product of the variance 10^{26*6}=10^{105}.
        Thus the result variance becomes overflow.
        '''
        adj = Adjugate(6, randRange=(1 << 17))
        ssVarOrg = addNoise(adj.ssOrg, 0)
        self.roundtrip(adj, 0, ssVarOrg)

        adj = Adjugate(6, randRange=(1 << 18))
        ssVarOrg = addNoise(adj.ssOrg, 0)
        with self.assertRaises((OverflowError, UncertaintyException)):
            self.roundtrip(adj, 0, ssVarOrg)


    def testConditionalNumber(self):
        '''
        detAdj = 10^{20}~10^{12}
        detRnd = 10^{102}~10^{176}
        '''
        adj = Adjugate(6, randRange=(1 << 12))
        ssVarOrg = addNoise(adj.ssOrg, 0)
        detAdj, detRnd, precStat = self.roundtrip(adj, 0, ssVarOrg)
        with self.assertRaises(LossUncertaintyException):
            detRnd / (detAdj**(5-1))


    @unittest.skip('Too slow')
    def testDump(self):
        '''
        Time to take with     
                MAX_SIZE = 10
                ELEMENT_RANGE = 1 << 9
                ROUND_RANGE = 1 << 9
            size=5, repeat=5, success: 0.083042
            size=6, repeat=4, success: 0.966857
            size=7, repeat=3, success: 12.42007
            size=8, repeat=0, fail: 95.956827
            size=4, repeat= 7, noise=1e-16: 0.195862
            size=5, repeat= 6, noise=1e-16: 1.473194
            size=6, repeat= 5, noise=1e-16: 19.911286
            size=7, repeat= 4, noise=1e-16: 255.788484
            size=8, repeat= 3, noise=1e-16: 3659.943234
        '''
        with open('./Python/Output/AdjMatrix.txt', 'w') as f, open('./Python/Output/DetMatrix.txt', 'w') as fc:
            f.write("NoiseType\tNoise\tSize\tRepeat\tLoss Uncertainty"
                    "\tError Deviation\tError Mean\tError Minimum\tError Maximum"
                    "\tForward Deviation\tForward Mean\tForward Minimum\tForward Maximum"
                    "\tRoundtrip Deviation\tRoundtrip Mean\tRoundtrip Minimum\tRoundtrip Maximum"
                    "\tUncertainty Deviation\tUncertainty Mean\tUncertainty Minimum\tUncertainty Maximum")
            for bucket in Adjugate.sAdjHisto[(Adjugate.MAX_SIZE - 1, 0)].buckets():
                f.write(f'\t{bucket:.1f}')
            f.write('\n')

            def write(adj, noise, repeat):
                f.write(f'Gaussian\t{noise}\t{adj.size}\t{repeat}')
                f.write(f'\t{Adjugate.sLossVar[(adj.size, noise)]/adj.size**2/2/(repeat + 1)}')
                histo = Adjugate.sAdjHisto[(adj.size, noise)]
                stat = histo.stat()
                f.write(f'\t{stat.dev()}\t{stat.mean()}\t{stat.min()}\t{stat.max()}')
                stat = Adjugate.sFwdStat[(adj.size, noise)]
                f.write(f'\t{stat.dev()}\t{stat.mean()}\t{stat.min()}\t{stat.max()}')
                stat = Adjugate.sRndStat[(adj.size, noise)]
                f.write(f'\t{stat.dev()}\t{stat.mean()}\t{stat.min()}\t{stat.max()}')
                stat = Adjugate.sAdjStat[(adj.size, noise)]
                f.write(f'\t{stat.dev()}\t{stat.mean()}\t{stat.min()}\t{stat.max()}')
                count = sum([c for c in histo.histogram()])
                for c in histo.histogram():
                        if count:
                            f.write(f'\t{c/count}')
                        else:
                            f.write(f'\t{c}')
                f.write('\n')
                f.flush()

            fc.write("NoiseType\tNoise\tSize\tRepeat"
                    "\tMatrix Value\tMatrix Uncertainty\tMatrix Precision"
                    "\tAdjugate Value\tAdjugate Uncertainty\tAdjugate Precision"
                    "\tPrecision Deviation\tPrecision Mean\tPrecision Minimum\tPrecision Maximum\n")

            def writeCond(adj, noise, repeat, detAdj, detRnd, stat):
                fc.write(f'Gaussian\t{noise}\t{adj.size}\t{repeat}')
                adjUnc = detAdj.uncertainty()
                fc.write(f'\t{detAdj.value()}\t{adjUnc}\t{adjUnc/abs(detAdj.value())}')
                rndUnc = detRnd.uncertainty()
                fc.write(f'\t{detRnd.value()}\t{rndUnc}\t{rndUnc/abs(detRnd.value())}')
                fc.write(f'\t{stat.dev()}\t{stat.mean()}\t{stat.min()}\t{stat.max()}\n')
                fc.flush()

            for size in range(5, Adjugate.MAX_SIZE):
                now = datetime.datetime.now()
                adj = Adjugate(size, randRange=(1 << 11))
                for repeat in range(Adjugate.MAX_SIZE - size):
                    ssVarOrg = addNoise(adj.ssOrg, 0)
                    try:
                        detAdj, detRnd, precStat = self.roundtrip(adj, 0, ssVarOrg)
                    except BaseException as ex:
                        try:
                            detAdj, detRnd, precStat = self.roundtrip(adj, 0, ssVarOrg)
                        except BaseException as ex:
                            print(f'Fail at size={size}: {ex}')
                            break
                    writeCond(adj, 0, repeat, detAdj, detRnd, precStat)
                else:
                    write(adj, 0, repeat)
                    print(f'size={size}, repeat={repeat}, success: {(datetime.datetime.now() - now).total_seconds()}')
                    continue
                print(f'size={size}, repeat={repeat}, fail: {(datetime.datetime.now() - now).total_seconds()}')
                break

            for size in range(4, Adjugate.MAX_SIZE):
                now = datetime.datetime.now()
                adj = Adjugate(size)
                for noise in NOISES:
                    for repeat in range(Adjugate.MAX_SIZE - size + 1):
                        ssVarOrg = addNoise(adj.ssOrg, Adjugate.noise(noise))
                        try:
                            detAdj, detRnd, precStat = self.roundtrip(adj, noise, ssVarOrg)
                        except BaseException as ex:
                            try:
                                detAdj, detRnd, precStat = self.roundtrip(adj, noise, ssVarOrg)
                            except BaseException as ex:
                                print(f'Fail at size={size} and noise={noise}: {ex}')
                                raise ex    
                        writeCond(adj, noise, repeat, detAdj, detRnd, precStat)
                    write(adj, noise, repeat)
                
                print(f'size={size}, repeat= {repeat}, noise={noise}: {(datetime.datetime.now() - now).total_seconds()}')

 


if __name__ == '__main__':
    unittest.main()