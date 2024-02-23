import math
import os
import unittest

from indexSin import IndexSin, RegressiveSin
from varDbl import VarDbl
from histo import Stat

class TestIndexSin (unittest.TestCase):

    def test_neg_rem(self):
        '''
        assume order==3, so size==8
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

    def test_get_sin_index(self):
        indexSin = IndexSin(3)
        self.assertEqual(8, indexSin.size())

        self.assertEqual( 0, indexSin.get_index(0))
        self.assertEqual( 1, indexSin.get_index(1))
        self.assertEqual( 2, indexSin.get_index(2))
        self.assertEqual( 3, indexSin.get_index(3))
        self.assertEqual( 4, indexSin.get_index(4))
        self.assertEqual( 3, indexSin.get_index(5))
        self.assertEqual( 2, indexSin.get_index(6))
        self.assertEqual( 1, indexSin.get_index(7))
        self.assertEqual( 0, indexSin.get_index(8))
        self.assertEqual(-1, indexSin.get_index(9))
        self.assertEqual(-2, indexSin.get_index(10))
        self.assertEqual(-3, indexSin.get_index(11))
        self.assertEqual(-4, indexSin.get_index(12))
        self.assertEqual(-3, indexSin.get_index(13))
        self.assertEqual(-2, indexSin.get_index(14))
        self.assertEqual(-1, indexSin.get_index(15))
        self.assertEqual( 0, indexSin.get_index(16))
        self.assertEqual( 1, indexSin.get_index(17))
        self.assertEqual( 2, indexSin.get_index(18))

        self.assertEqual(-1, indexSin.get_index(-1))
        self.assertEqual(-2, indexSin.get_index(-2))
        self.assertEqual(-3, indexSin.get_index(-3))
        self.assertEqual(-4, indexSin.get_index(-4))
        self.assertEqual(-3, indexSin.get_index(-5))
        self.assertEqual(-2, indexSin.get_index(-6))
        self.assertEqual(-1, indexSin.get_index(-7))
        self.assertEqual( 0, indexSin.get_index(-8))
        self.assertEqual( 1, indexSin.get_index(-9))
        self.assertEqual( 2, indexSin.get_index(-10))
        self.assertEqual( 3, indexSin.get_index(-11))
        self.assertEqual( 4, indexSin.get_index(-12))
        self.assertEqual( 3, indexSin.get_index(-13))
        self.assertEqual( 2, indexSin.get_index(-14))
        self.assertEqual( 1, indexSin.get_index(-15))
        self.assertEqual( 0, indexSin.get_index(-16))
        self.assertEqual(-1, indexSin.get_index(-17))
        self.assertEqual(-2, indexSin.get_index(-18))

    def test_sin(self):
        indexSin = IndexSin(3)
        self.assertEqual(8, indexSin.size())

        self.assertEqual(indexSin.sin(0), 0)
        self.assertEqual(indexSin.sin(1), math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(2), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(3), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(4), 1)
        self.assertEqual(indexSin.sin(5), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(6), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(7), math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(8), 0)

        self.assertEqual(indexSin.sin(-1), -math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(-2), -math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(-3), -math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(-4), -1)
        self.assertEqual(indexSin.sin(-5), -math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(-6), -math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(-7), -math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(-8), 0)

    def test_sin_error(self):
        self.assertEqual(-math.sin(1/8*math.pi) - 9.992007221626409e-16, 
                         math.sin(1001/8*math.pi))
        self.assertEqual(math.ulp(math.sin(1/8*math.pi)), 5.551115123125783e-17)

    def test_cos(self):
        indexSin = IndexSin(3)
        self.assertEqual(8, indexSin.size())

        self.assertEqual(indexSin.cos(0), 1)
        self.assertEqual(indexSin.cos(1), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.cos(2), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.cos(3), math.sin(1/8*math.pi))
        self.assertEqual(indexSin.cos(4), 0)
        self.assertEqual(indexSin.cos(5), -math.sin(1/8*math.pi))
        self.assertEqual(indexSin.cos(6), -math.cos(2/8*math.pi))
        self.assertEqual(indexSin.cos(7), -math.cos(1/8*math.pi))
        self.assertEqual(indexSin.cos(8), -1)

        self.assertEqual(indexSin.cos(-1), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.cos(-2), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.cos(-3), math.sin(1/8*math.pi))

    def test_arc_sin(self):
        indexSin = IndexSin(3)
        self.assertEqual( 0, indexSin.arc_sin( 0))
        self.assertEqual( 1, indexSin.arc_sin( math.sin(math.pi*1/8)))
        self.assertEqual( 2, indexSin.arc_sin( math.sin(math.pi*2/8)))
        self.assertEqual( 3, indexSin.arc_sin( math.sin(math.pi*3/8)))
        self.assertEqual( 4, indexSin.arc_sin( 1))

        self.assertEqual(-1, indexSin.arc_sin(-math.sin(math.pi*1/8)))
        self.assertEqual(-2, indexSin.arc_sin(-math.sin(math.pi*2/8)))
        self.assertEqual(-3, indexSin.arc_sin(-math.sin(math.pi*3/8)))
        self.assertEqual(-4, indexSin.arc_sin(-1))

        self.assertEqual( 2.239966192595619, indexSin.arc_sin( math.sin(math.pi*1/3)))
        self.assertEqual(-2.239966192595619, indexSin.arc_sin(-math.sin(math.pi*1/3)))


class TestIndexSinVsLibSin (unittest.TestCase):

    @unittest.skip('Too slow')
    def testSin(self):
        with open(f'./Python/Output/SinDiff.txt', 'w') as f:
            f.write('Order\tFreq\tX\tIndex Sin\tValue Error\tUncertainty\n')
            for order in range(4, 8):
                size = 1 << order
                sin = IndexSin(order)
                for i in range(size**2 //2):
                    x = i / size * math.pi
                    err = VarDbl(math.sin(x)) - VarDbl(sin.sin(i))
                    f.write(f'{order}\t{i}\t{x}\t{sin.sin(i)}\t{err.value()}\t{err.uncertainty()}\n')  

        with open(f'./Python/Output/SinDiff_Stat.txt', 'w') as f:
            f.write('Order\tPart\tCount\tMin\tMax\tAbsMax\tMean\tDev\n')
            for order in range(4, 17):
                sin = IndexSin(order)
                size = 1 << order

                stat = Stat()
                for i in range(size >> 3):
                    x = i / size * math.pi
                    err = VarDbl(math.sin(x))**2 + VarDbl(math.cos(x))**2 - 1
                    stat.accum( err.value() / err.uncertainty() )
                absMax = max(abs(stat.min()), abs(stat.max()))
                f.write(f'{order}\tIndexed\t{stat.count()}\t{stat.min()}\t{stat.max()}\t{absMax}\t{stat.mean()}\t{stat.dev()}\n')

                stat = Stat()
                for i in range(size >> 3, size**2 //2):
                    x = i / size * math.pi
                    stat.accum( math.sin(x) - sin.sin(i) )
                absMax = max(abs(stat.min()), abs(stat.max()))
                f.write(f'{order}\tLib\t{stat.count()}\t{stat.min()}\t{stat.max()}\t{absMax}\t{stat.mean()}\t{stat.dev()}\n')

                f.flush()
       

    @unittest.skip('Too slow')
    def testCtan(self):
        with open(f'./Python/Output/CotDiff.txt', 'w') as f:
            f.write('Order\tFreq\tX\tIndex Cot\tValue Error\tUncertainty\n')
            for order in range(4, 8):
                size = 1 << order
                sin = IndexSin(order)
                for i in range(size):
                    x = i / size * math.pi
                    try:
                        ctan = VarDbl(sin.cos(i)) / VarDbl(sin.sin(i))
                        err = ctan - 1/math.tan(x)
                        f.write(f'{order}\t{i}\t{x}\t{ctan.value()}\t{err.value()}\t{err.uncertainty()}\n') 
                    except BaseException as ex:
                        print(f'Ignore i={i} order={order}: {ex}')
                        continue

            for order in range(8, 17):
                size = 1 << order
                sin = IndexSin(order)
                sin = IndexSin(order)
                for i in range(size *7//8, size):
                    x = i / size * math.pi
                    try:
                        ctan = VarDbl(sin.cos(i)) / VarDbl(sin.sin(i))
                        err = ctan - 1/math.tan(x)
                        f.write(f'{order}\t{i}\t{x}\t{ctan.value()}\t{err.uncertainty()}\t{err.value()}\n') 
                    except BaseException as ex:
                        print(f'Ignore i={i} order={order}: {ex}')
                        continue




class TestRegressiveSin(unittest.TestCase):

    def testWithUncertainty_4(self):
        sin = RegressiveSin(4)
        sin.calc()
        self.assertIsNone(sin.withUncertainty())

    def testWithUncertainty_6(self):
        sin = RegressiveSin(6)
        sin.calc()
        self.assertIsNone(sin.withUncertainty())

    def testWithUncertainty_18(self):
        sin = RegressiveSin(18)
        if not os.path.isfile(RegressiveSin.path(18)):
            sin.calc()
        self.assertIsNone(sin.withUncertainty())

    @unittest.skip("too slow")
    def testWithUncertainty_19(self):
        sin = RegressiveSin(19)
        if not os.path.isfile(RegressiveSin.path(19)):
            sin.calc()
        self.assertIsNone(sin.withUncertainty())


if __name__ == '__main__':
    unittest.main()