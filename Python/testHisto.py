import math
import unittest

from histo import Stat, Histo

class TestStat (unittest.TestCase):

    def assertEmpty(self, stat:Stat):
        self.assertEqual(0, stat.count())
        self.assertFalse(math.isfinite(stat.min()))
        self.assertFalse(math.isfinite(stat.max()))
        self.assertFalse(math.isfinite(stat.mean()))
        self.assertFalse(math.isfinite(stat.dev()))

    def testInit(self):
        stat = Stat()
        self.assertEmpty(stat)

    def testAccumFailed(self):
        stat = Stat()
        stat.accum(float('nan'))
        self.assertEmpty(stat)
        stat.accum(float('inf'))
        self.assertEmpty(stat)
        stat.accum(float('-inf'))
        self.assertEmpty(stat)

    def testAccum(self):
        stat = Stat()
        self.assertTrue(stat.accum(0))
        self.assertEqual(1, stat.count())
        self.assertEqual(0, stat.min())
        self.assertEqual(0, stat.max())
        self.assertEqual(0, stat.mean())
        self.assertEqual(0, stat.dev())   # not stddev
        self.assertEqual('"Stat: 1, 0.0"', str(stat))

        self.assertTrue(stat.accum(1))
        self.assertTrue(stat.accum(-1))
        self.assertEqual(3, stat.count())
        self.assertEqual(-1, stat.min())
        self.assertEqual(+1, stat.max())
        self.assertEqual(0, stat.mean())
        self.assertEqual(math.sqrt(2/3), stat.dev())
        self.assertEqual('"Stat: 3, 0.0+/-0.816496580927726"', str(stat))


class TestHisto (unittest.TestCase):

    def testRange(self):
        histo = Histo(2, 1.5)
        self.assertEqual(1.5, histo.range())
        histo = Histo(2, 1.6)
        self.assertEqual(1.5, histo.range())
        histo = Histo(2, 1.4)
        self.assertEqual(1.5, histo.range())

    def testHisto(self):
        histo = Histo(2, 1.5)
        self.assertListEqual([-3/2,-1,-1/2,0.0,1/2,1,3/2], histo.buckets())
        for i in range(-4,5):
            self.assertTrue(histo.accum(i/2))
        self.assertEqual(9, histo.stat().count())
        self.assertEqual(-2, histo.stat().min())
        self.assertEqual(+2, histo.stat().max())
        self.assertEqual(0,  histo.stat().mean())
        self.assertEqual(math.sqrt(5/3),  histo.stat().dev())
        self.assertEqual(1,  histo.less())
        self.assertEqual(1,  histo.more())
        self.assertListEqual([1]*7, histo.histogram())
        self.assertEqual('"Stat: 9, 0.0+/-1.2909944487358056"', str(histo))

    def testAdd(self):
        histo = Histo(2, 1.5)
        self.assertListEqual([-3/2,-1,-1/2,0.0,1/2,1,3/2], histo.buckets())
        for i in range(-4,0):
            self.assertTrue(histo.accum(i/2))
        histo1 = Histo(2, 1.5)
        self.assertListEqual([-3/2,-1,-1/2,0.0,1/2,1,3/2], histo1.buckets())
        for i in range(0,5):
            self.assertTrue(histo1.accum(i/2))
        histo += histo1
        self.assertEqual(9, histo.stat().count())
        self.assertEqual(-2, histo.stat().min())
        self.assertEqual(+2, histo.stat().max())
        self.assertEqual(0,  histo.stat().mean())
        self.assertEqual(math.sqrt(5/3),  histo.stat().dev())
        self.assertEqual(1,  histo.less())
        self.assertEqual(1,  histo.more())
        self.assertListEqual([1]*7, histo.histogram())
        self.assertEqual('"Stat: 9, 0.0+/-1.2909944487358056"', str(histo))




if __name__ == '__main__':
    unittest.main()