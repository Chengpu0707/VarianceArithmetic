import datetime
from fractions import Fraction
import itertools
import logging
import math
import numpy
import os
import unittest

from histo import Histo, Stat
from matrix import permutSign, isSquareMatrix, createIntMatrix, createHilbertMatrix, addNoise
from matrix import linear, multiply, adjugate, adjugate_mul
from taylor import NotReliableException, NotMonotonicException, NotFiniteException
from varDbl import VarDbl, InitException


logger = logging.getLogger(__name__)



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

    def testNumpyIntMatrix(self):
        '''
        Can not use numpy to hold int matrix and to retain int for determinant
        '''
        matrix = createIntMatrix(2, randRange=Adjugate.ELEMENT_RANGE)
        self.assertTrue(isSquareMatrix(matrix, sType=(int,)))
        mat = numpy.matrix(matrix)
        self.assertEqual(mat.dtype, numpy.int64)
        self.assertEqual(type(numpy.linalg.det(mat)), numpy.float64)

    def testCreateIntMatrix(self):
        '''
        Test deviation of the random int matrix 
        '''
        valStat = Stat()
        for size in range(4, Adjugate.MAX_SIZE):
            ssOrg = createIntMatrix(size, randRange=Adjugate.ELEMENT_RANGE)
            self.assertTrue(isSquareMatrix(ssOrg, sType=(int,)))
            try:
                for i in range(size):
                    for j in range(size):
                        valStat.accum(ssOrg[i][j])
                dev =  Adjugate.noise(1)
                self.assertTrue(abs(valStat.mean()) < dev/2, f'size {size}: val_mean = {valStat.mean()} > {dev/2}, val_dev = {valStat.dev()}')
                self.assertAlmostEqual(valStat.dev(), dev, delta=dev/4, msg=f'size {size}: val_mean = {valStat.mean()}, val_dev = {valStat.dev()} vs {dev}')
            except BaseException as ex:
                raise ex

    def testAddNoise(self):
        '''
        Verify the input noise
        '''
        noise = 1e-3
        for size in range(4, Adjugate.MAX_SIZE):
            dev = Adjugate.noise(noise)
            ssOrg = createIntMatrix(size, randRange=Adjugate.ELEMENT_RANGE)
            ssMat = addNoise(ssOrg, dev)
            self.assertTrue(isSquareMatrix(ssMat, sType=(VarDbl,)))
            uncStat = Stat()
            errStat = Stat()
            nrmStat = Stat()
            try:
                for i in range(size):
                    for j in range(size):
                        self.assertEqual(type(ssMat[i][j]), VarDbl)
                        self.assertEqual(type(ssMat[i][j].value()), float)
                        unc = ssMat[i][j].uncertainty()
                        uncStat.accum(unc)
                        errStat.accum(ssMat[i][j].value() - ssOrg[i][j])
                        nrmStat.accum((ssMat[i][j].value() - ssOrg[i][j])/unc)
                self.assertEqual(uncStat.dev(), 0)
                self.assertAlmostEqual(uncStat.mean(), dev, 
                                       msg=f'size {size}: uncertainty_mean = {uncStat.mean()}, uncertainty_dev = {uncStat.dev()} vs {dev}')
                self.assertAlmostEqual(errStat.dev(), dev, delta=dev/2, 
                                       msg=f'size {size}: error_mean = {errStat.mean()} > {dev/2}, error_dev = {errStat.dev()}')
                self.assertTrue(0.5 <= nrmStat.dev() <= 1.5, 
                                msg=f'size {size}: normalized_dev = {nrmStat.dev()} not in [0.5, 1.5]')
                self.assertTrue(-0.5 <= nrmStat.mean() <= 0.5, 
                                msg=f'size {size}: normalized_mean = {nrmStat.mean()} not in [-0.5, +0.5]')
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


# The following constants belongs to Adjugate, but they can not be defined in Adjugate
NOISES = tuple([0] + [math.pow(10, nexp) for nexp in range(-17, 0)])


class Adjugate (unittest.TestCase):
    '''
    A class to hold precise matrix for adjugate matrix.
    A precise matrix has all int elements.
    ELEMENT_RANGE decides MAX_SIZE so that the calculated variance does not overflow
    '''
    ELEMENT_RANGE = 1 << 8
    MIN_SIZE = 2
    MAX_SIZE = 9    # overflow at matrix size 9 for ELEMENT_RANGE = 1 << 8

    @staticmethod
    def noise(noiseLevel:float) ->float:
        return Adjugate.ELEMENT_RANGE/math.sqrt(3) * noiseLevel

    sAdjHisto = {(size, noise): Histo(5, 3) 
                 for size in range(MAX_SIZE) for noise in NOISES}
    sAdjStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in NOISES}
    sAdjLoss = {(size, noise): 0 for size in range(MAX_SIZE) for noise in NOISES}
    sFwdStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in NOISES}
    sFwdLoss = {(size, noise): 0 for size in range(MAX_SIZE) for noise in NOISES}
    sRndStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in NOISES}
    sRndLoss = {(size, noise): 0 for size in range(MAX_SIZE) for noise in NOISES}
    sMulStat = {(size, noise): Stat() for size in range(MAX_SIZE) for noise in NOISES}
    sMulLoss = {(size, noise): 0 for size in range(MAX_SIZE) for noise in NOISES}

    __slots__ = ('size', 'randRange', 'ssOrg', 'ssAdj', 'detAdj')

    def __init__(self, size:int, randRange=ELEMENT_RANGE) -> None:
        '''
        Create a random int matrix of given size {self.ssOrg}
        Calculate its adjugate int matrix {self.ssAdj} and its determinant {self.detAdj}
        '''
        super().__init__()

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
                    self.assertEqual(ssId[i][j], (self.detAdj if i == j else 0))
                    self.assertEqual(ssId2[i][j], (detRnd if i == j else 0))
                    self.assertEqual(ssRnd[i][j] * self.detAdj, self.ssOrg[i][j] * detRnd)
                except AssertionError as ex:
                    raise ex
                
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

    def verifyIdentity(self, det, ssId, delta=1e-6, uncRange=3, places=5):
        size = len(ssId)
        for i in range(size):
            for j in range(size):
                try:
                    eType = type(ssId[i][j])
                    if eType == VarDbl:
                        self.assertAlmostEqual(ssId[i][j].value(), det.value() if i == j else 0,
                                               delta=max(abs(det.value())*delta, ssId[i][j].uncertainty()*uncRange))
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

    def roundtrip(self, noise:float, 
                  verify:bool=False, countDeterminant:bool=False):
        '''
        Carry out roundtrip test for given adjugate matrix {adj} and input matrix {ssOrg}
        The input matrix ssOrg is expected to be close to self.ssOrg with noise level noise

        "verify"=True may fail for matrix size larger than 6
        "countDeterminant"=True may change the stats for Adjugate.sRndStat to depends on matrix size.
        '''
        ssOrg = addNoise(self.ssOrg, Adjugate.noise(noise))
        detAdj, ssAdj = adjugate(ssOrg)
        ssId = multiply(ssOrg, ssAdj)      
        try:
            self.verifyIdentity(detAdj, multiply(ssOrg, ssAdj))
        except AssertionError as ex:
            if verify:
                raise ex
            else:
                logger.warning(f'Adjugate det={detAdj}, ex={ex}: {ssOrg}')

        detRnd, ssRnd = adjugate(ssAdj)
        try:
            self.verifyIdentity(detRnd, multiply(ssAdj, ssRnd), places=5)
        except AssertionError as ex:
            if verify:
                raise ex
            else:
                logger.warning(f'Roundtrip det={detRnd}, ex={ex}: {ssOrg}')

        detMul, ssMul = adjugate_mul(ssOrg)
        try:
            self.assertAlmostEqual(detMul.value(), detAdj.value())
        except AssertionError as ex:
            if verify:
                raise ex
            else:
                logger.warning(f'Multiply detMul={detMul} vs detAdj={detAdj} ex={ex}: {ssOrg}')

        def accum(value, expected, adjStat:Stat, adjHisto:Histo):
            if type(value) == VarDbl:
                adjUnc = value.uncertainty()
                adjStat.accum(adjUnc)
                if adjUnc:
                    adjHisto.accum((value.value() - expected)/adjUnc)
            else:
                Adjugate.sAdjLoss[(self.size, noise)] += 1

        if countDeterminant:
            adjStat = Adjugate.sAdjStat[(self.size + 1, noise)]
            adjHisto = Adjugate.sAdjHisto[(self.size + 1, noise)]
            accum(detAdj, self.detAdj, adjStat, adjHisto)

        adjStat = Adjugate.sAdjStat[(self.size, noise)]
        adjHisto = Adjugate.sAdjHisto[(self.size, noise)]

        for i in range(self.size):
            for j in range(self.size):
                accum(ssAdj[i][j], self.ssAdj[i][j], adjStat, adjHisto)

                diff = ssId[i][j] - detAdj if i == j else ssId[i][j]
                if (type(diff) == VarDbl) and diff.uncertainty():
                    Adjugate.sFwdStat[(self.size, noise)].accum(diff.value()/diff.uncertainty())
                else:
                    Adjugate.sFwdLoss[(self.size, noise)] += 1

                diff = detAdj * ssRnd[i][j] - detRnd * ssOrg[i][j]
                if (type(diff) == VarDbl) and diff.uncertainty():
                    Adjugate.sRndStat[(self.size, noise)].accum(diff.value()/diff.uncertainty())
                else:
                    Adjugate.sRndLoss[(self.size, noise)] += 1

                try:
                    self.assertAlmostEqual(ssMul[i][j].value(), ssAdj[i][j].value() if type(ssAdj[i][j]) == VarDbl else ssAdj[i][j])
                except AssertionError as ex:
                    if verify:
                        raise ex
                    else:
                        logger.warning(f'Multiply ssMul[{i}][{j}]={ssMul[i][j]} vs ssAdj[{i}][{j}]={ssAdj[i][j]} ex={ex}: {ssOrg}')
                diff = ssMul[i][j] - self.ssAdj[i][j]
                if (type(diff) == VarDbl) and diff.uncertainty():
                    Adjugate.sMulStat[(self.size, noise)].accum(diff.value()/diff.uncertainty())
                else:
                    Adjugate.sMulLoss[(self.size, noise)] += 1

        return detAdj
    
    @staticmethod
    def dumpPath(minSize:int=MIN_SIZE, maxSize:int=MAX_SIZE) ->tuple[str, str]:
        if os.getcwd().endswith('\\VarianceArithmetic'):
            dumpPath = f'./Python/Output/AdjMatrix_{minSize}_{maxSize}.txt'
            logPath = f'./Python/Output/AdjMatrix_{minSize}_{maxSize}.log'
        elif os.getcwd().endswith('\\VarianceArithmetic\\Python'):
            dumpPath = f'./Output/AdjMatrix_{minSize}_{maxSize}.txt'
            logPath = f'./Output/AdjMatrix_{minSize}_{maxSize}.log'
        else:
            raise RuntimeError(f'Wrong working directory: {os.getcwd()}')
        return dumpPath, logPath
    
    @staticmethod
    def dump(minSize:int=MIN_SIZE, maxSize:int=MAX_SIZE, sampleCount:int=1024):
        dumpPath, logPath = Adjugate.dumpPath(minSize, maxSize)
        if os.path.isfile(logPath):
            try:
                os.remove(logPath)
            except:
                pass
        logging.basicConfig(filename=logPath, encoding='utf-8', level=logging.DEBUG,
                            format='%(asctime)s:%(levelname)s:%(message)s')
        
        def  statHeader(name:str) ->str:
            return f'\t{name} ' + f'\t{name} '.join(('Deviation', 'Mean', 'Minimum', 'Maximum', 'Count', 'Loss'))
        
        def writeStat(stat:Stat, loss:int):
            f.write(f'\t{stat.dev()}\t{stat.mean()}\t{stat.min()}\t{stat.max()}\t{stat.count()}\t{loss}')

        header = "NoiseType\tNoise\tSize\tRepeat" + \
                 statHeader('Adjugate Error') + statHeader('Adjugate Uncertainty') + \
                 statHeader('Forward Error') + statHeader('Roundtrip Error') + \
                 statHeader('Multiple Error') + '\t' + \
                 '\t'.join(str(i) for i in Adjugate.sAdjHisto[(maxSize - 1, 0)].buckets()) + '\n' 
        columnCount = header.count('\t') + 1       
        
        sExist = {}
        exist = os.path.isfile(dumpPath)
        if exist:
            with open(dumpPath) as f:
                hdr = next(f)
                if hdr != header:
                    raise RuntimeError(f'Invalid header {hdr} vs {header}')
                for ln, line in enumerate(f):
                    sWord = line.split('\t')
                    if (columnCount != len(sWord)) or (sWord[0] != 'Gaussian'):
                        raise RuntimeError(f'Invalid line {ln + 1}: {line}')
                    noise = float(sWord[1])
                    if noise not in NOISES:
                        raise RuntimeError(f'Invalid line {ln + 1} noise {noise} not in {NOISES}: {line}')
                    size = int(sWord[2])
                    if size < minSize or size >= maxSize:
                        raise RuntimeError(f'Invalid line {ln + 1} size {size} not in [{minSize}, {maxSize}): {line}')
                    repeat = int(sWord[3])
                    if repeat < 0:
                        raise RuntimeError(f'Invalid line {ln + 1} repeat {repeat}: {line}')
                    adjDev = float(sWord[4])
                    if adjDev < 0:
                        raise RuntimeError(f'Invalid line {ln + 1} adjDev {adjDev}: {line}')
                    sExist[(size, noise)] = adjDev

        
        with open(dumpPath, 'a' if exist else 'w') as f:
            if not exist:
                f.write(f'{header}\n')

            def write(adj, noise, repeat):
                f.write(f'Gaussian\t{noise}\t{adj.size}\t{repeat}')
                histo = Adjugate.sAdjHisto[(adj.size, noise)]
                writeStat(histo.stat(), histo.less() + histo.more())
                writeStat(Adjugate.sAdjStat[(adj.size, noise)], Adjugate.sAdjLoss[(adj.size, noise)])
                writeStat(Adjugate.sFwdStat[(adj.size, noise)], Adjugate.sFwdLoss[(adj.size, noise)])
                writeStat(Adjugate.sRndStat[(adj.size, noise)], Adjugate.sRndLoss[(adj.size, noise)])
                writeStat(Adjugate.sMulStat[(adj.size, noise)], Adjugate.sMulLoss[(adj.size, noise)])
                count = sum([c for c in histo.histogram()])
                if count:
                    for c in histo.histogram():
                        f.write(f'\t{c/count}')
                else:
                    f.write('\t'.join([''] * len(histo.buckets())))
                f.write('\n')
                f.flush()

            for size in range(minSize, maxSize):
                for noise in NOISES:
                    if (size, noise) in sExist:
                        logger.info(f'Skip size={size}, noise={noise} already exists with adjDev={sExist[(size, noise)]}')
                        continue
                    print(f'{datetime.datetime.now()}: Start size={size}, noise={noise}')
                    for repeat in range(sampleCount // size**2):
                        adj = Adjugate(size)
                        logger.info(f'Start  size={size}, noise={noise}, repeat={repeat}, detAdj={adj.detAdj}: {adj.ssOrg}')
                        try:
                            detAdj = adj.roundtrip(noise, ssVarOrg)
                        except BaseException as ex:  # avoid singular
                            logger.info(f'First failure to process size={size}, noise={noise}, repeat={repeat}: ex={ex}')
                            try:
                                adj = Adjugate(size)
                                logger.info(f'Start size={size}, noise={noise}, repeat={repeat}, detAdj={adj.detAdj}: {adj.ssOrg}')
                                ssVarOrg = addNoise(adj.ssOrg, Adjugate.noise(noise))
                                detAdj = adj.roundtrip(noise)
                            except BaseException as ex:
                                logger.warning(f'Second failure to process size={size}, noise={noise}, repeat={repeat}: ex={ex}')
                                raise ex
                        logger.info(f'Finish size={size}, noise={noise}, repeat={repeat}, detAdj={detAdj}: {ssVarOrg}')
                    write(adj, noise, repeat)
    



class TestAdjugate (unittest.TestCase):

    def verifyValue(self, val, ret):
        Adjugate.verifyValue(self, val, ret)

    def verifyIdentity(self, det, ssId, delta=1e-6, uncRange=3, places=5):
        Adjugate.verifyIdentity(self, det, ssId, delta, uncRange, places)

    def verify(self, ssMat, det, ssAdj):
        Adjugate.verify(self, ssMat, det, ssAdj)

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
        ssMat = ((VarDbl(1,1e-1), VarDbl(2, 1e-2)), 
                 (VarDbl(3,1e-3), VarDbl(4, 1e-4)))
        det = VarDbl(1,1e-1)*VarDbl(4,1e-4) - VarDbl(2,1e-2)*VarDbl(3,1e-3)
        ssAdj = ((VarDbl( 4,1e-4), VarDbl(-2,1e-2)), 
                 (VarDbl(-3,1e-3), VarDbl( 1,1e-1)))
        self.verify(ssMat, det, ssAdj)


    def testIntSize3(self):
        self.verify(((1,0,0), (0,1,0), (0,0,1)), 1, ((1,0,0), (0,1,0), (0,0,1)))
        self.verify(((1,2,3), (-4,-5,6), (7,8,9)), 72, ((-93,6,27), (78,-12,-18), (3,6,3)))

    def testVarDblSize3_noTracing(self):
        noTracing = VarDbl( 1,1e-1)*VarDbl(-5,1e-5)*VarDbl( 9,1e-9) - VarDbl( 1,1e-1)*VarDbl( 6,1e-6)*VarDbl( 8,1e-8) +\
                    VarDbl( 2,1e-2)*VarDbl( 6,1e-6)*VarDbl( 7,1e-7) - VarDbl( 2,1e-2)*VarDbl(-4,1e-4)*VarDbl( 9,1e-9) +\
                    VarDbl( 3,1e-3)*VarDbl(-4,1e-4)*VarDbl( 8,1e-8) - VarDbl( 3,1e-3)*VarDbl(-5,1e-5)*VarDbl( 7,1e-7)
        self.verifyValue(noTracing, VarDbl(72, 6.602897701208691))

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
                VarDbl(-93, math.sqrt((1e-5*1e-9)**2 + (5*1e-9)**2 + (1e-5*9)**2 +
                                      (1e-6*1e-8)**2 + (6*1e-8)**2 + (1e-6*8)**2)),
                VarDbl(6,   math.sqrt((1e-2*1e-9)**2 + (2*1e-9)**2 + (1e-2*9)**2 +
                                      (1e-8*1e-3)**2 + (8*1e-3)**2 + (1e-8*3)**2)),
                VarDbl(27,  math.sqrt((1e-2*1e-6)**2 + (2*1e-5)**2 + (1e-2*6)**2 +
                                      (1e-5*1e-3)**2 + (5*1e-3)**2 + (1e-5*3)**2)),
            ), 
            (
                VarDbl(78,  math.sqrt((1e-5*1e-7)**2 + (6*1e-7)**2 + (1e-5*7)**2 +
                                      (1e-4*1e-9)**2 + (4*1e-9)**2 + (1e-4*9)**2)),
                VarDbl(-12, math.sqrt((1e-1*1e-9)**2 + (1*1e-9)**2 + (1e-1*9)**2 +
                                      (1e-3*1e-7)**2 + (3*1e-7)**2 + (1e-3*7)**2)),
                VarDbl(-18, math.sqrt((1e-1*1e-6)**2 + (1*1e-6)**2 + (1e-1*6)**2 +
                                      (1e-4*1e-3)**2 + (4*1e-3)**2 + (1e-4*3)**2)),
            ),
            (
                VarDbl(3,   math.sqrt((1e-4*1e-8)**2 + (4*1e-8)**2 + (1e-4*8)**2 +
                                      (1e-5*1e-7)**2 + (5*1e-7)**2 + (1e-5*7)**2)),
                VarDbl(6,   math.sqrt((1e-2*1e-7)**2 + (2*1e-7)**2 + (1e-2*7)**2 + 
                                      (1e-1*1e-8)**2 + (1*1e-8)**2 + (1e-1*8)**2)),
                VarDbl(3,   math.sqrt((1e-1*1e-5)**2 + (1*1e-5)**2 + (1e-1*5)**2 +
                                      (1e-2*1e-4)**2 + (2*1e-4)**2 + (1e-2*4)**2)),
            ))
        self.verify(ssMat3, VarDbl(val3, math.sqrt(var3)), ssAdj)


    def testDirectMultiplication(self):
        ssMat3 = ((VarDbl( 1,1e-1), VarDbl( 2,1e-2), VarDbl( 3,1e-3)),
                  (VarDbl(-4,1e-4), VarDbl(-5,1e-5), VarDbl( 6,1e-6)),
                  (VarDbl( 7,1e-7), VarDbl( 8,1e-8), VarDbl( 9,1e-9)))
        det, ssMul = adjugate_mul(ssMat3)
        self.assertAlmostEqual(det.value(), 72)
        self.assertListEqual([ssMul[i][j] for i in range(3) for j in range(3)],
                             [-93,6,27, 78, -12,-18, 3,6,3] )


    def testException(self):
        '''
        The adjugate matrix has order of magnitude of 2^{18*6} == 10^{31}.  
        The actual measured is 10^{28}~10^{26}.
        The roundtrip matrix has order of magnitude of 10^{28*6} == 10^{206}.  
        The square of the value is 10^{412}. The product of the variance 10^{26*6}=10^{105}.
        Thus the result variance becomes overflow.
        This test may occationally fail due to how random matrix is generated
        '''
        adj = Adjugate(6, randRange=(1 << 17))
        adj.roundtrip(0)

        adj = Adjugate(6, randRange=(1 << 19))
        with self.assertRaises((OverflowError, InitException)):
            adj.roundtrip(0.1)


    def testHilbert(self):
        '''
        numpy can not accept float
        '''
        ssIntHilbert = createHilbertMatrix(6)
        with self.assertRaises(BaseException):
            numpy.linalg.cond(numpy.asarray(ssIntHilbert))
        cond = numpy.linalg.cond(numpy.asarray([[float(val) for val in row] for row in ssIntHilbert]))
        self.assertAlmostEqual(cond, 14951058.6424659)
        ssHilbert = addNoise(ssIntHilbert, 0)
        det, ssAdj = adjugate(ssHilbert)
        self.assertAlmostEqual(numpy.linalg.cond(numpy.asarray([[val.value() for val in row] for row in ssHilbert])), 
                               cond)
        self.verifyIdentity(det, multiply(ssHilbert, ssAdj), delta=1e-1, uncRange=6)


    def testInversionException(self):
        ssOrgVal = (
            ( 154.6037896556856,    -241.41721719977542,    -215.5313557001085), 
            (  62.317618798899076,  -173.31255384513335,     255.8635167560795), 
            (-198.35031352268214,    225.851083348709,       141.43828027978913)
        )
        ssOrg = tuple([tuple([VarDbl(val, 14.780166891254423) for val in sOrgVal]) for sOrgVal in ssOrgVal])
        detAdj, ssAdj = adjugate(ssOrg)
        self.assertAlmostEqual(detAdj.value(),       6031779.536578279)
        self.assertAlmostEqual(detAdj.uncertainty(), 2311975.810766203)
        try:
            ssAdj[0][0] / detAdj
            self.fail()
        except NotMonotonicException:
            pass
        except NotReliableException:
            pass
        
        adj = Adjugate(2)
        for noise in (0, 1e-2, 1e-1):
            ssOrg = addNoise(adj.ssOrg, Adjugate.noise(noise))
            detAdj, ssAdj = adjugate(ssOrg)
            for i in range(adj.size):
                for j in range(adj.size):
                    try:
                        (ssAdj[i][j] /detAdj) - (adj.ssAdj[i][j] /adj.detAdj)
                    except (NotMonotonicException, NotReliableException, NotFiniteException) as ex:
                        print(f'Found at ssAdj[{i}][{j}]={ssAdj[i][j]}, detAdj={detAdj}, noise={noise}, ssOrg={ssOrg}')
                        break
        

    def testInversionValue(self):
        '''
        Calculation for (((1, 1e-1), (2, 1e-2)), ((3, 1e-3), (4, 1e-4)))^-1
        '''
        common = ((1**2)*((1e-4)**2) + (4**2)*((1e-1)**2) + (2**2)*((1e-3)**2) + (3**2)*((1e-2)**2))/(2**4)
        self.assertAlmostEqual(common, 0.010056500625000003)

        sX = [[4/-2 - 1/(2**2)*(1e-4)**2, -2/-2 - -3/(2**2)*(1e-2)**2], 
              [-3/-2 - -2/(2**2)*(1e-3)**2, 1/-2 - 4/(2**2)*(1e-1)**2]]
        self.assertAlmostEqual(sX[0][0], -2.0000000025)
        self.assertAlmostEqual(sX[0][1], 1.000075)
        self.assertAlmostEqual(sX[1][0], 1.5000005)
        self.assertAlmostEqual(sX[1][1], -0.51)

        sdX = [
            [ math.sqrt(((1e-4)**2)*((1/(2**2))+4*1*4/(2**3))+(4**2)*common),
              math.sqrt(((1e-2)**2)*((1/(2**2))+4*2*3/(2**3))+(2**2)*common) ],
            [ math.sqrt(((1e-3)**2)*((1/(2**2))+4*2*3/(2**3))+(3**2)*common),
              math.sqrt(((1e-1)**2)*((1/(2**2))+4*1*4/(2**3))+(1**2)*common) ]
        ]
        self.assertAlmostEqual(sdX[0][0], 0.40112844887890964)
        self.assertAlmostEqual(sdX[0][1], 0.2013727948358467)
        self.assertAlmostEqual(sdX[1][0], 0.30085171700523833)
        self.assertAlmostEqual(sdX[1][1], 0.1804342002642515)

        dir = VarDbl(1, 0.1) / VarDbl(-2, 0.161)
        self.assertAlmostEqual(dir.value(), -0.5033051985276559)
        self.assertAlmostEqual(dir.uncertainty(), 0.06526322323589419)


    def testAdjugate(self):
        Adjugate.dump(2, 6)
        dumpPath, logPath = Adjugate.dumpPath(2, 6)
        self.assertTrue(os.path.isfile(dumpPath))
        os.remove(dumpPath)
        self.assertTrue(os.path.isfile(logPath))


    def testIdealCoverage(self):
        '''
        May fail occassionally due to statstical stablility
        '''
        for size in range(4, 7):
            adj = Adjugate(size)
            for noise in (1e-3, 1e-9):
                for repeat in range(Adjugate.MAX_SIZE - size):
                    adj.roundtrip(noise, verify=True)

                adjHisto = Adjugate.sAdjHisto[(adj.size, noise)]
                rndStat = Adjugate.sRndStat[(adj.size, noise)]
                try:
                    self.assertGreater(adjHisto.stat().dev(), 0.5)
                    self.assertLess(adjHisto.stat().dev(), 2)
                    self.assertLess(abs(adjHisto.stat().mean()), 0.3)
                    self.assertLess(rndStat.dev(), 0.1)
                    self.assertLess(rndStat.mean(), 0.1)
                    self.assertGreater(Adjugate.sMulStat[(adj.size, noise)].dev(), 0.5)
                    self.assertLess(Adjugate.sMulStat[(adj.size, noise)].dev(), 2)
                except BaseException as ex:
                    raise ex

    def testProperCoverage(self):
        for size in range(4, 7):
            adj = Adjugate(size, randRange=(1 << 14))
            for repeat in range(Adjugate.MAX_SIZE - size):
                adj.roundtrip(0)
            adjHisto = Adjugate.sAdjHisto[(adj.size, 0)]
            if adjHisto.stat().dev():
                try:
                    self.assertGreater(adjHisto.stat().dev(), 0.1)
                    self.assertLess(adjHisto.stat().mean(), 10)
                except BaseException as ex:
                    raise ex


                


if __name__ == '__main__':
    unittest.main()