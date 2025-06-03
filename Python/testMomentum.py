import math
import os
import unittest

import momentum


class TestNormal (unittest.TestCase):

    def testNormal(self):
        mmt = momentum.Normal()
        self.assertEqual(mmt.maxOrder, 442)
        self.assertAlmostEqual(mmt[0],    1*(1 - 5.733031e-07))
        self.assertEqual(mmt[1], 0)
        self.assertAlmostEqual(mmt[2],    1*(1 - 1.544050e-05))
        self.assertEqual(mmt[3], 0)
        self.assertAlmostEqual(mmt[4],    3*(1 - 1.393338e-04))
        self.assertEqual(mmt[5], 0)
        self.assertAlmostEqual(mmt[6],   15*(1 - 7.588003e-04))
        self.assertEqual(mmt[7], 0)
        self.assertAlmostEqual(mmt[8],  105*(1 - 2.9711805e-03))
        self.assertEqual(mmt[9], 0)
        self.assertAlmostEqual(mmt[10], 945*(1 - 9.1166811e-03))

    def testCalcLow(self):
        self.assertAlmostEqual(momentum.Normal.calcLow(5, 0),    1*(1 - 5.733031e-07))
        self.assertEqual(momentum.Normal.calcLow(5, 1), 0)
        self.assertAlmostEqual(momentum.Normal.calcLow(5, 2),    1*(1 - 1.544050e-05))
        self.assertEqual(momentum.Normal.calcLow(5, 3), 0)
        self.assertAlmostEqual(momentum.Normal.calcLow(5, 4),    3*(1 - 1.393338e-04))
        self.assertEqual(momentum.Normal.calcLow(5, 5), 0)
        self.assertAlmostEqual(momentum.Normal.calcLow(5, 6),   15*(1 - 7.588003e-04))
        self.assertEqual(momentum.Normal.calcLow(5, 7), 0)
        self.assertAlmostEqual(momentum.Normal.calcLow(5, 8),  105*(1 - 2.9711805e-03))
        self.assertEqual(momentum.Normal.calcLow(5, 9), 0)
        self.assertAlmostEqual(momentum.Normal.calcLow(5, 10), 945*(1 - 9.1166811e-03))

    def testApproxCalc(self):
        mmt = momentum.Normal()
        self.assertTrue(math.isfinite(mmt._sMomentum[-1]))
        mmt2 = momentum.Normal(divid=mmt._divid*2, readCached=False)
        self.assertTrue(math.isfinite(mmt2._sMomentum[-1]))
        self.assertAlmostEqual(mmt2._sMomentum[0] / mmt._sMomentum[0], 1)
        self.assertAlmostEqual(mmt2._sMomentum[1] / mmt._sMomentum[1], 1)
        self.assertAlmostEqual(mmt2._sMomentum[2] / mmt._sMomentum[2], 1)
        self.assertAlmostEqual(mmt2._sMomentum[3] / mmt._sMomentum[3], 1)
        self.assertAlmostEqual(mmt2._sMomentum[4] / mmt._sMomentum[4], 1)
        self.assertAlmostEqual(mmt2._sMomentum[5] / mmt._sMomentum[5], 1, delta=6e-8)
        self.assertAlmostEqual(mmt2._sMomentum[6] / mmt._sMomentum[6], 1, delta=2e-7)
        self.assertAlmostEqual(mmt2._sMomentum[7] / mmt._sMomentum[7], 1, delta=2e-7)
        self.assertAlmostEqual(mmt2._sMomentum[8] / mmt._sMomentum[8], 1, delta=3e-7)
        self.assertAlmostEqual(mmt2._sMomentum[9] / mmt._sMomentum[9], 1, delta=4e-7)
        self.assertAlmostEqual(mmt2._sMomentum[10] / mmt._sMomentum[10], 1, delta=4e-7)
        for i in range(20, len(mmt._sMomentum)):
            self.assertAlmostEqual(mmt2._sMomentum[i] / mmt._sMomentum[i], 1, delta=5e-5)

    def testMaxOrder(self):
        mmt = momentum.Normal(readCached=False)
        self.assertEqual(mmt.maxOrder, 442)
        self.assertEqual(len(mmt._sMomentum), mmt.maxOrder // 2)
        self.assertAlmostEqual(mmt._sMomentum[-1], 1.2365430968160818e+300)


if __name__ == '__main__':
    unittest.main()