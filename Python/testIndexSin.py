import math
import os
import typing
import unittest

from indexSin import IndexSin, SinSource, OUTDIR
from varDbl import VarDbl
from histo import Stat, Histo



q1 = math.sin(1/8*math.pi)
q2 = math.sin(2/8*math.pi)
q3 = math.cos(1/8*math.pi)


class TestException (unittest.TestCase):

    def testSize(self):
        self.assertEqual(IndexSin.validateSize(1024), 10)
        with self.assertRaises(RuntimeError):
            IndexSin.validateSize(1000)

    def testorder(self):
        for order in range(IndexSin.MIN_ORDER, IndexSin.MAX_ORDER + 1):
            IndexSin.validateOrder(order)
        with self.assertRaises(RuntimeError):
            IndexSin.validateSize(-1)
        with self.assertRaises(RuntimeError):
            IndexSin.validateSize(IndexSin.MAX_ORDER + 1)


class TestIndexSin (unittest.TestCase):

    def test_float_precision(self):
        self.assertEqual(float('0.70710678118654757274'), 0.7071067811865476)
        self.assertEqual(float('7.07106781186547572737e-01'), 0.7071067811865476)

    def test_neg_rem(self):
        '''
        The reminder is always positive
        '''
        # whe index is near 0
        self.assertEqual(0, 1 // 4)     #0 % 1
        self.assertEqual(1, 1 % 4)
        self.assertEqual(-1, -1 // 4)   #0 % -1
        self.assertEqual(3, -1 % 4)

        # whe index is near +pi/2
        self.assertEqual(0, 3 // 4)     #0 % 3
        self.assertEqual(3, 3 % 4)
        self.assertEqual(1, 5 // 4)     #0 % 3
        self.assertEqual(1, 5 % 4)

        # whe index is near +pi
        self.assertEqual(1, 7 // 4)     #0 % 1
        self.assertEqual(3, 7 % 4)
        self.assertEqual(2, 9 // 4)     #0 % -1
        self.assertEqual(1, 9 % 4)

        # whe index is near -pi/2
        self.assertEqual(-1, -3 // 4)   #0 % -3
        self.assertEqual(1, -3 % 4)
        self.assertEqual(-2, -5 // 4)   #0 % -3
        self.assertEqual(3, -5 % 4)

        # whe index is near -pi
        self.assertEqual(-2, -7 // 4)   #0 % -1
        self.assertEqual(1, -7 % 4)
        self.assertEqual(-3, -9 // 4)   #0 % 1
        self.assertEqual(3, -9 % 4)


    def assert_get_index_1(self, indexSin):
        self.assertEqual( 0, indexSin.get_index(0, 1))
        self.assertEqual( 1, indexSin.get_index(1, 1))
        self.assertEqual( 0, indexSin.get_index(2, 1))
        self.assertEqual(-1, indexSin.get_index(3, 1))
        self.assertEqual( 0, indexSin.get_index(4, 1))

        self.assertEqual(-1, indexSin.get_index(-1, 1))
        self.assertEqual( 0, indexSin.get_index(-2, 1))
        self.assertEqual( 1, indexSin.get_index(-3, 1))
        self.assertEqual( 0, indexSin.get_index(-4, 1))

    def assert_get_index_2(self, indexSin):
        self.assertEqual( 0, indexSin.get_index(0, 2))
        self.assertEqual( 1, indexSin.get_index(1, 2))
        self.assertEqual( 2, indexSin.get_index(2, 2))
        self.assertEqual( 3, indexSin.get_index(3, 2))
        self.assertEqual( 0, indexSin.get_index(4, 2))
        self.assertEqual(-1, indexSin.get_index(5, 2))
        self.assertEqual(-2, indexSin.get_index(6, 2))
        self.assertEqual(-3, indexSin.get_index(7, 2))
        self.assertEqual( 0, indexSin.get_index(8, 2))
        self.assertEqual( 1, indexSin.get_index(9, 2))

        self.assertEqual(-3, indexSin.get_index(-1, 2))
        self.assertEqual(-2, indexSin.get_index(-2, 2))
        self.assertEqual(-1, indexSin.get_index(-3, 2))
        self.assertEqual( 0, indexSin.get_index(-4, 2))
        self.assertEqual( 3, indexSin.get_index(-5, 2))
        self.assertEqual( 2, indexSin.get_index(-6, 2))
        self.assertEqual( 1, indexSin.get_index(-7, 2))
        self.assertEqual( 0, indexSin.get_index(-8, 2))
        self.assertEqual(-3, indexSin.get_index(-9, 2))

    def assert_get_index_3(self, indexSin):
        self.assertEqual( 0, indexSin.get_index(0, 3))
        self.assertEqual( 1, indexSin.get_index(1, 3))
        self.assertEqual( 2, indexSin.get_index(2, 3))
        self.assertEqual( 3, indexSin.get_index(3, 3))
        self.assertEqual( 4, indexSin.get_index(4, 3))
        self.assertEqual( 5, indexSin.get_index(5, 3))
        self.assertEqual( 6, indexSin.get_index(6, 3))
        self.assertEqual( 7, indexSin.get_index(7, 3))
        self.assertEqual( 0, indexSin.get_index(8, 3))
        self.assertEqual(-1, indexSin.get_index(9, 3))
        self.assertEqual(-2, indexSin.get_index(10, 3))
        self.assertEqual(-3, indexSin.get_index(11, 3))
        self.assertEqual(-4, indexSin.get_index(12, 3))
        self.assertEqual(-5, indexSin.get_index(13, 3))
        self.assertEqual(-6, indexSin.get_index(14, 3))
        self.assertEqual(-7, indexSin.get_index(15, 3))
        self.assertEqual( 0, indexSin.get_index(16, 3))
        self.assertEqual( 1, indexSin.get_index(17, 3))
        self.assertEqual( 2, indexSin.get_index(18, 3))

        self.assertEqual(-7, indexSin.get_index(-1, 3))
        self.assertEqual(-6, indexSin.get_index(-2, 3))
        self.assertEqual(-5, indexSin.get_index(-3, 3))
        self.assertEqual(-4, indexSin.get_index(-4, 3))
        self.assertEqual(-3, indexSin.get_index(-5, 3))
        self.assertEqual(-2, indexSin.get_index(-6, 3))
        self.assertEqual(-1, indexSin.get_index(-7, 3))
        self.assertEqual( 0, indexSin.get_index(-8, 3))
        self.assertEqual( 7, indexSin.get_index(-9, 3))
        self.assertEqual( 6, indexSin.get_index(-10, 3))
        self.assertEqual( 5, indexSin.get_index(-11, 3))
        self.assertEqual( 4, indexSin.get_index(-12, 3))
        self.assertEqual( 3, indexSin.get_index(-13, 3))
        self.assertEqual( 2, indexSin.get_index(-14, 3))
        self.assertEqual( 1, indexSin.get_index(-15, 3))
        self.assertEqual( 0, indexSin.get_index(-16, 3))
        self.assertEqual(-7, indexSin.get_index(-17, 3))
        self.assertEqual(-6, indexSin.get_index(-18, 3))


    def assert_get_half_index_1(self, indexSin):
        self.assertEqual( 0, indexSin.get_index(0, 1))
        self.assertEqual( 1, indexSin.get_index(1, 1))
        self.assertEqual( 0, indexSin.get_index(2, 1))
        self.assertEqual(-1, indexSin.get_index(3, 1))
        self.assertEqual( 0, indexSin.get_index(4, 1))

        self.assertEqual(-1, indexSin.get_index(-1, 1))
        self.assertEqual( 0, indexSin.get_index(-2, 1))
        self.assertEqual( 1, indexSin.get_index(-3, 1))
        self.assertEqual( 0, indexSin.get_index(-4, 1))

    def assert_get_half_index_2(self, indexSin):
        self.assertEqual( 0, indexSin.get_index(0, 2))
        self.assertEqual( 1, indexSin.get_index(1, 2))
        self.assertEqual( 2, indexSin.get_index(2, 2))
        self.assertEqual( 1, indexSin.get_index(3, 2))
        self.assertEqual( 0, indexSin.get_index(4, 2))
        self.assertEqual(-1, indexSin.get_index(5, 2))
        self.assertEqual(-2, indexSin.get_index(6, 2))
        self.assertEqual(-1, indexSin.get_index(7, 2))
        self.assertEqual( 0, indexSin.get_index(8, 2))
        self.assertEqual( 1, indexSin.get_index(9, 2))

        self.assertEqual(-1, indexSin.get_index(-1, 2))
        self.assertEqual(-2, indexSin.get_index(-2, 2))
        self.assertEqual(-1, indexSin.get_index(-3, 2))
        self.assertEqual( 0, indexSin.get_index(-4, 2))
        self.assertEqual( 1, indexSin.get_index(-5, 2))
        self.assertEqual( 2, indexSin.get_index(-6, 2))
        self.assertEqual( 1, indexSin.get_index(-7, 2))
        self.assertEqual( 0, indexSin.get_index(-8, 2))
        self.assertEqual(-1, indexSin.get_index(-9, 2))

    def assert_get_half_index_3(self, indexSin):
        self.assertEqual( 0, indexSin.get_index(0, 3))
        self.assertEqual( 1, indexSin.get_index(1, 3))
        self.assertEqual( 2, indexSin.get_index(2, 3))
        self.assertEqual( 3, indexSin.get_index(3, 3))
        self.assertEqual( 4, indexSin.get_index(4, 3))
        self.assertEqual( 3, indexSin.get_index(5, 3))
        self.assertEqual( 2, indexSin.get_index(6, 3))
        self.assertEqual( 1, indexSin.get_index(7, 3))
        self.assertEqual( 0, indexSin.get_index(8, 3))
        self.assertEqual(-1, indexSin.get_index(9, 3))
        self.assertEqual(-2, indexSin.get_index(10, 3))
        self.assertEqual(-3, indexSin.get_index(11, 3))
        self.assertEqual(-4, indexSin.get_index(12, 3))
        self.assertEqual(-3, indexSin.get_index(13, 3))
        self.assertEqual(-2, indexSin.get_index(14, 3))
        self.assertEqual(-1, indexSin.get_index(15, 3))
        self.assertEqual( 0, indexSin.get_index(16, 3))
        self.assertEqual( 1, indexSin.get_index(17, 3))
        self.assertEqual( 2, indexSin.get_index(18, 3))

        self.assertEqual(-1, indexSin.get_index(-1, 3))
        self.assertEqual(-2, indexSin.get_index(-2, 3))
        self.assertEqual(-3, indexSin.get_index(-3, 3))
        self.assertEqual(-4, indexSin.get_index(-4, 3))
        self.assertEqual(-3, indexSin.get_index(-5, 3))
        self.assertEqual(-2, indexSin.get_index(-6, 3))
        self.assertEqual(-1, indexSin.get_index(-7, 3))
        self.assertEqual( 0, indexSin.get_index(-8, 3))
        self.assertEqual( 1, indexSin.get_index(-9, 3))
        self.assertEqual( 2, indexSin.get_index(-10, 3))
        self.assertEqual( 3, indexSin.get_index(-11, 3))
        self.assertEqual( 4, indexSin.get_index(-12, 3))
        self.assertEqual( 3, indexSin.get_index(-13, 3))
        self.assertEqual( 2, indexSin.get_index(-14, 3))
        self.assertEqual( 1, indexSin.get_index(-15, 3))
        self.assertEqual( 0, indexSin.get_index(-16, 3))
        self.assertEqual(-1, indexSin.get_index(-17, 3))
        self.assertEqual(-2, indexSin.get_index(-18, 3))


    def assert_sin_1(self, indexSin):
        self.assertAlmostEqual(indexSin.sin(0, 1).value(), 0)
        self.assertAlmostEqual(indexSin.sin(1, 1).value(), 1)
        self.assertAlmostEqual(indexSin.sin(2, 1).value(), 0)
        self.assertAlmostEqual(indexSin.sin(3, 1).value(), -1)
        self.assertAlmostEqual(indexSin.sin(4, 1).value(), 0)

        self.assertAlmostEqual(indexSin.sin(-1, 1).value(), -1)
        self.assertAlmostEqual(indexSin.sin(-2, 1).value(), 0)

    def assert_sin_2(self, indexSin):
        self.assertAlmostEqual(indexSin.sin(0, 2).value(), 0)
        self.assertAlmostEqual(indexSin.sin(1, 2).value(), q2)
        self.assertAlmostEqual(indexSin.sin(2, 2).value(), 1)
        self.assertAlmostEqual(indexSin.sin(3, 2).value(), q2)
        self.assertAlmostEqual(indexSin.sin(4, 2).value(), 0)
        self.assertAlmostEqual(indexSin.sin(5, 2).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(6, 2).value(), -1)
        self.assertAlmostEqual(indexSin.sin(7, 2).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(8, 2).value(), 0)

        self.assertAlmostEqual(indexSin.sin(-1, 2).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(-2, 2).value(), -1)
        self.assertAlmostEqual(indexSin.sin(-3, 2).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(-4, 2).value(), 0)

    def assert_sin_3(self, indexSin):
        self.assertAlmostEqual(indexSin.sin(0, 3).value(), 0)
        self.assertAlmostEqual(indexSin.sin(1, 3).value(), q1)
        self.assertAlmostEqual(indexSin.sin(2, 3).value(), q2)
        self.assertAlmostEqual(indexSin.sin(3, 3).value(), q3)
        self.assertAlmostEqual(indexSin.sin(4, 3).value(), 1)
        self.assertAlmostEqual(indexSin.sin(5, 3).value(), q3)
        self.assertAlmostEqual(indexSin.sin(6, 3).value(), q2)
        self.assertAlmostEqual(indexSin.sin(7, 3).value(), q1)
        self.assertAlmostEqual(indexSin.sin(8, 3).value(), 0)
        self.assertAlmostEqual(indexSin.sin(9, 3).value(), -q1)
        self.assertAlmostEqual(indexSin.sin(10, 3).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(11, 3).value(), -q3)
        self.assertAlmostEqual(indexSin.sin(12, 3).value(), -1)
        self.assertAlmostEqual(indexSin.sin(13, 3).value(), -q3)
        self.assertAlmostEqual(indexSin.sin(14, 3).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(15, 3).value(), -q1)
        self.assertAlmostEqual(indexSin.sin(16, 3).value(), 0)

        self.assertAlmostEqual(indexSin.sin(-1, 3).value(), -q1)
        self.assertAlmostEqual(indexSin.sin(-2, 3).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(-3, 3).value(), -q3)
        self.assertAlmostEqual(indexSin.sin(-4, 3).value(), -1)
        self.assertAlmostEqual(indexSin.sin(-5, 3).value(), -q3)
        self.assertAlmostEqual(indexSin.sin(-6, 3).value(), -q2)
        self.assertAlmostEqual(indexSin.sin(-7, 3).value(), -q1)
        self.assertAlmostEqual(indexSin.sin(-8, 3).value(), 0)


    def assert_cos_1(self, indexSin):
        self.assertAlmostEqual(indexSin.cos(0, 1).value(), 1)
        self.assertAlmostEqual(indexSin.cos(1, 1).value(), 0)
        self.assertAlmostEqual(indexSin.cos(2, 1).value(), -1)
        self.assertAlmostEqual(indexSin.cos(3, 1).value(), 0)
        self.assertAlmostEqual(indexSin.cos(4, 1).value(), 1)

        self.assertAlmostEqual(indexSin.cos(-1, 1).value(), 0)

    def assert_cos_2(self, indexSin):
        self.assertAlmostEqual(indexSin.cos(0, 2).value(), 1)
        self.assertAlmostEqual(indexSin.cos(1, 2).value(), q2)
        self.assertAlmostEqual(indexSin.cos(2, 2).value(), 0)
        self.assertAlmostEqual(indexSin.cos(3, 2).value(), -q2)
        self.assertAlmostEqual(indexSin.cos(4, 2).value(), -1)
        self.assertAlmostEqual(indexSin.cos(5, 2).value(), -q2)
        self.assertAlmostEqual(indexSin.cos(6, 2).value(), 0)
        self.assertAlmostEqual(indexSin.cos(7, 2).value(), q2)
        self.assertAlmostEqual(indexSin.cos(8, 2).value(), 1)

        self.assertAlmostEqual(indexSin.cos(-1, 2).value(), q2)
        self.assertAlmostEqual(indexSin.cos(-2, 2).value(), 0)

    def assert_cos_3(self, indexSin):
        self.assertAlmostEqual(indexSin.cos(0, 3).value(), 1)
        self.assertAlmostEqual(indexSin.cos(1, 3).value(), q3)
        self.assertAlmostEqual(indexSin.cos(2, 3).value(), q2)
        self.assertAlmostEqual(indexSin.cos(3, 3).value(), q1)
        self.assertAlmostEqual(indexSin.cos(4, 3).value(), 0)
        self.assertAlmostEqual(indexSin.cos(5, 3).value(), -q1)
        self.assertAlmostEqual(indexSin.cos(6, 3).value(), -q2)
        self.assertAlmostEqual(indexSin.cos(7, 3).value(), -q3)
        self.assertAlmostEqual(indexSin.cos(8, 3).value(), -1)
        self.assertAlmostEqual(indexSin.cos(9, 3).value(), -q3)
        self.assertAlmostEqual(indexSin.cos(10, 3).value(), -q2)
        self.assertAlmostEqual(indexSin.cos(11, 3).value(), -q1)
        self.assertAlmostEqual(indexSin.cos(12, 3).value(), 0)
        self.assertAlmostEqual(indexSin.cos(13, 3).value(), q1)
        self.assertAlmostEqual(indexSin.cos(14, 3).value(), q2)
        self.assertAlmostEqual(indexSin.cos(15, 3).value(), q3)
        self.assertAlmostEqual(indexSin.cos(16, 3).value(), 1)

        self.assertAlmostEqual(indexSin.cos(-1, 3).value(), q3)
        self.assertAlmostEqual(indexSin.cos(-2, 3).value(), q2)
        self.assertAlmostEqual(indexSin.cos(-3, 3).value(), q1)


    def writeToFile(self, indexSin:IndexSin, order:int):
        indexSin.dump(order)


class TestPrec (TestIndexSin):
    indexSin = IndexSin(SinSource.Prec)

    def test_index(self):
        self.assert_get_half_index_1(TestPrec.indexSin)
        self.assert_get_half_index_2(TestPrec.indexSin)
        self.assert_get_half_index_3(TestPrec.indexSin)

    def test_sin(self):
        self.assert_sin_1(TestPrec.indexSin)
        self.assert_sin_2(TestPrec.indexSin)
        self.assert_sin_3(TestPrec.indexSin)

    def test_cos(self):
        self.assert_cos_1(TestPrec.indexSin)
        self.assert_cos_2(TestPrec.indexSin)
        self.assert_cos_3(TestPrec.indexSin)


class TestPrecAdj (TestIndexSin):
    indexSin = IndexSin(SinSource.PrecAdj)

    def test_index(self):
        self.assert_get_half_index_1(TestPrecAdj.indexSin)
        self.assert_get_half_index_2(TestPrecAdj.indexSin)
        self.assert_get_half_index_3(TestPrecAdj.indexSin)

    def test_sin(self):
        self.assert_sin_1(TestPrecAdj.indexSin)
        self.assert_sin_2(TestPrecAdj.indexSin)
        self.assert_sin_3(TestPrecAdj.indexSin)

    def test_cos(self):
        self.assert_cos_1(TestPrecAdj.indexSin)
        self.assert_cos_2(TestPrecAdj.indexSin)
        self.assert_cos_3(TestPrecAdj.indexSin)

    def test_writeToFile(self):
        self.writeToFile(TestPrecAdj.indexSin, 10)

    def test_writeToFile_large(self):
        self.writeToFile(TestPrecAdj.indexSin, 18)


class TestQuart (TestIndexSin):
    indexSin = IndexSin()

    def test_index(self):
        self.assert_get_half_index_1(TestQuart.indexSin)
        self.assert_get_half_index_2(TestQuart.indexSin)
        self.assert_get_half_index_3(TestQuart.indexSin)

    def test_sin(self):
        self.assert_sin_1(TestQuart.indexSin)
        self.assert_sin_2(TestQuart.indexSin)
        self.assert_sin_3(TestQuart.indexSin)

    def test_cos(self):
        self.assert_cos_1(TestQuart.indexSin)
        self.assert_cos_2(TestQuart.indexSin)
        self.assert_cos_3(TestQuart.indexSin)

    def test_writeToFile(self):
        self.writeToFile(TestQuart.indexSin, 10)

    def test_writeToFile_large(self):
        self.writeToFile(TestQuart.indexSin, 18)



class TestFull (TestIndexSin):
    indexSin = IndexSin(SinSource.Full)

    def test_get_index(self):
        self.assert_get_index_1(TestFull.indexSin)
        self.assert_get_index_2(TestFull.indexSin)
        self.assert_get_index_3(TestFull.indexSin)

    def test_sin(self):
        self.assert_sin_1(TestFull.indexSin)
        self.assert_sin_2(TestFull.indexSin)
        self.assert_sin_3(TestFull.indexSin)

    def test_cos(self):
        self.assert_cos_3(TestFull.indexSin)

    def test_writeToFile(self):
        self.writeToFile(TestFull.indexSin, 10)

    def test_writeToFile_large(self):
        self.writeToFile(TestFull.indexSin, 18)


class TestFixed (TestIndexSin):
    indexSin = IndexSin(SinSource.Fixed)

    def test_get_index(self):
        self.assert_get_index_1(TestFixed.indexSin)
        self.assert_get_index_2(TestFixed.indexSin)
        self.assert_get_index_3(TestFixed.indexSin) 

    def test_sin(self):
        self.assert_sin_1(TestFixed.indexSin)
        self.assert_sin_2(TestFixed.indexSin)
        self.assert_sin_3(TestFixed.indexSin)

    def test_cos(self):
        self.assert_cos_3(TestFixed.indexSin)

    def test_writeToFile_large(self):
        self.writeToFile(TestFixed.indexSin, 18)


class TestLimit (TestIndexSin):
    indexSin = IndexSin(SinSource.Limit, 
                        sCosSin=[VarDbl(v) for v in (1,0, q3,q1, q2,q2, q1,q3, 0,1, -q1,q3, -q2,q2, -q3,q1)])

    def test_get_index(self):
        self.assert_get_index_1(TestLimit.indexSin) 
        self.assert_get_index_2(TestLimit.indexSin) 
        self.assert_get_index_3(TestLimit.indexSin) 

    def test_sin(self):
        self.assert_sin_1(TestLimit.indexSin)
        self.assert_sin_2(TestLimit.indexSin)
        self.assert_sin_3(TestLimit.indexSin)

    def test_cos(self):
        TestIndexSin.assert_cos_3(self, TestLimit.indexSin)

class TestLib (TestIndexSin):
    indexSin = IndexSin(SinSource.Lib)

    def test_sin(self):
        self.assert_sin_1(TestLib.indexSin)
        self.assert_sin_2(TestLib.indexSin)
        self.assert_sin_3(TestLib.indexSin)

    def test_cos(self):
        TestIndexSin.assert_cos_3(self, TestLib.indexSin)

    def test_writeToFile_large(self):
        self.writeToFile(TestLib.indexSin, 18)


def dumpStat(dumpPath:str, testName:str, maxRange:float,
             sStat:dict[int, Stat], sHist:dict[int, Histo]):
    '''
    dump the uncertainty stat {sStat} and normalize error histogram {sHist} to file {dumpPath}
    '''
    if (not sHist) or set(sStat.keys()) != set(sHist.keys()):
        raise RuntimeError(f'Invalid uncertainty stat {sStat.keys()} and normalize error histogram {sHist.keys()}')
    with open(dumpPath, 'w') as fw:
        hist = sHist[list(sHist.keys())[0]]
        HEADER = ("Order\tRange\tTest"
                  "\tUncertainty Count\tUncertainty Mean\tUncertainty Deviation\tUncertainty Minimum\tUncertainty Maximum"
                  "\tError Count\tError Mean\tError Deviation\tError Minimum\tError Maximum\t") + \
                  '\t'.join([str(i) for i in hist.buckets()])
        fw.write(HEADER)
        for order in sorted(sHist):
            stat = sStat[order]
            hist = sHist[order]
            fw.write(f'\n{order}\t{maxRange}\t{testName}\t'
                        f'{stat.count()}\t{stat.mean()}\t{stat.dev()}\t{stat.min()}\t{stat.max()}\t'
                        f'{hist.stat().count()}\t{hist.stat().mean()}\t{hist.stat().dev()}\t{hist.stat().min()}\t{hist.stat().max()}\t')
            if hist.stat().count():
                fw.write('\t'.join([str(c/hist.stat().count()) for c in hist.histogram()]))
            else:
                fw.write('\t' *(len(hist.buckets()) - 1))
            fw.flush()


class TestIndexSinDiff (unittest.TestCase):
    '''
    Compare between different SinSource of IndexSin
    '''
    quart = IndexSin(SinSource.Quart)
    full = IndexSin(SinSource.Full)
    prec = IndexSin(SinSource.Prec)
    precAdj = IndexSin(SinSource.PrecAdj)

    def profile(self, testName:str, 
                funcIdx:typing.Callable[[int, int], VarDbl], funcLib:typing.Callable[[int, int], VarDbl],
                sOrderRange=(4,11), maxRange=8,
                divids=5, devs=3):
        '''
        Compare two sin(i, order) between [0, {maxRange}] for the order in {sOrder}
        The output is display with x = i/(1 << order), without PI 
        '''
        sStat = {}
        sHist = {}
        dumpPath = f'{OUTDIR}/Python/Output/Profile_{testName}_{sOrderRange[0]}_{sOrderRange[-1]}_{maxRange}.txt'
        HEADER = "Order\tRange\tTest\tX\tValue Error\tUncertainty\tNormalized Error"
        with open(dumpPath, 'w') as fw:
            fw.write(HEADER)
            for order in range(*sOrderRange):
                sStat[order] = Stat()
                sHist[order] = Histo(divids=divids, devs=devs)
                size = 1 << order
                for i in range(size*maxRange + 1):
                    x = i / size
                    fw.write(f'\n{order}\t{maxRange}\t{testName}\t{x}')
                    try:
                        valIdx = funcIdx(i, order)
                        valLib = funcLib(i, order)
                        err = valLib - valIdx
                        sStat[order].accum(err.uncertainty())
                        if err.uncertainty() > 0:
                            sHist[order].accum(err.value() / err.uncertainty())
                        fw.write(f'\t{err.value()}\t{err.uncertainty()}')
                        if err.uncertainty() > 0:
                            fw.write(f'\t{err.value() /err.uncertainty()}')
                        else:
                            fw.write('\t')
                    except BaseException as ex:
#                        traceback.print_exc()
                        continue
                    fw.flush()

        dumpPath = f'{OUTDIR}/Python/Output/Stat_{testName}_{sOrderRange[0]}_{sOrderRange[-1]}_{maxRange}.txt'
        dumpStat(dumpPath, testName, maxRange, sStat, sHist)

    def test_SinDiff_Quart_vs_Prec(self):
        self.profile('SinDiff_Quart_vs_Prec',
                     lambda i, order: TestIndexSinDiff.quart.sin(i, order),
                     lambda i, order: TestIndexSinDiff.prec.sin(i, order))
        
    def test_SinDiff_Prec_vs_Lib(self):
        self.profile('SinDiff_Prec_vs_Lib', 
                     lambda i, order: TestIndexSinDiff.prec.sin(i, order),
                     lambda i, order: VarDbl(math.sin(math.pi *i/(1 << order))))

    def test_SinDiff_Quart_vs_Lib(self):
        self.profile('SinDiff_Quart_vs_Lib',
                     lambda i, order: TestIndexSinDiff.quart.sin(i, order),
                     lambda i, order: VarDbl(math.sin(math.pi *i/(1 << order))))
        
    def test_SinDiff_Quart_vs_Full(self):
        self.profile('SinDiff_Quart_vs_Full', 
                     lambda i, order: TestIndexSinDiff.quart.sin(i, order),
                     lambda i, order: TestIndexSinDiff.full.sin(i, order),
                     maxRange=1)


    def test_SinError_Quart(self):
        self.profile('SinError_Quart', 
                     lambda i, order: TestIndexSinDiff.quart.sin(i, order)**2 + TestIndexSinDiff.quart.cos(i, order)**2,
                     lambda i, order: 1)
        
    def test_SinError_Lib(self):
        self.profile('SinError_Lib', 
                     lambda i, order: 1,
                     lambda i, order: VarDbl(math.sin(math.pi *i/(1 << order)))**2 + VarDbl(math.cos(math.pi *i/(1 << order)))**2)

    def test_SinError_Prec(self):
        self.profile('SinError_Prec', 
                     lambda i, order: TestIndexSinDiff.prec.sin(i, order)**2 + TestIndexSinDiff.prec.cos(i, order)**2,
                     lambda i, order: 1)

    def test_SinError_PrecAdj(self):
        self.profile('SinError_PrecAdj', 
                     lambda i, order: TestIndexSinDiff.precAdj.sin(i, order)**2 + TestIndexSinDiff.precAdj.cos(i, order)**2,
                     lambda i, order: 1)


    def test_CTanDiff_Quart_vs_Lib(self):
        self.profile('CTanDiff_Quart_vs_Lib', 
                     lambda i, order: TestIndexSinDiff.quart.cos(i, order) / TestIndexSinDiff.quart.sin(i, order),
                     lambda i, order: VarDbl(1 /math.tan(math.pi *i/(1 << order))))
        
    def test_CTanDiff_Lib_vs_Lib(self):
        self.profile('CTanDiff_Lib_vs_Lib', 
                     lambda i, order: VarDbl(math.cos(math.pi *i/(1 << order)) / math.sin(math.pi *i/(1 << order))),
                     lambda i, order: VarDbl(1 /math.tan(math.pi *i/(1 << order))))

    def test_CTanDiff_Quart_vs_Full(self):
        self.profile('CTanDiff_Quart_vs_Full', 
                     lambda i, order: TestIndexSinDiff.quart.cos(i, order) / TestIndexSinDiff.quart.sin(i, order),
                     lambda i, order: TestIndexSinDiff.full.cos(i, order) / TestIndexSinDiff.full.sin(i, order),
                     maxRange=1)
        

    def test_TanDiff_Quart_vs_Lib(self):
        self.profile('TanDiff_Quart_vs_Lib', 
                     lambda i, order: TestIndexSinDiff.quart.sin(i, order) / TestIndexSinDiff.quart.cos(i, order),
                     lambda i, order: VarDbl(math.tan(math.pi *i/(1 << order))))
        
    def test_TanDiff_Lib_vs_Lib(self):
        self.profile('TanDiff_Lib_vs_Lib', 
                     lambda i, order: VarDbl(math.sin(math.pi *i/(1 << order)) / math.cos(math.pi *i/(1 << order))),
                     lambda i, order: VarDbl(math.tan(math.pi *i/(1 << order))))
        
    def test_TanDiff_Quart_vs_Full(self):
        self.profile('TanDiff_Quart_vs_Full', 
                     lambda i, order: TestIndexSinDiff.quart.sin(i, order) / TestIndexSinDiff.quart.cos(i, order),
                     lambda i, order: TestIndexSinDiff.full.sin(i, order) / TestIndexSinDiff.full.cos(i, order),
                     maxRange=1)








if __name__ == '__main__':
    unittest.main()