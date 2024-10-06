import math
import os
import unittest

import momentum


class TestMomentum (unittest.TestCase):

    def testNormal(self):
        mmt = momentum.Momentum()
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

    def testApproxCalc(self):
        if os.path.isfile(momentum.Momentum.FILE_APPROX[0]):
            os.remove(momentum.Momentum.FILE_APPROX[0])
        elif os.path.isfile(momentum.Momentum.FILE_APPROX[1]):
            os.remove(momentum.Momentum.FILE_APPROX[1])
        with self.assertRaises(ValueError):
            mmt = momentum.Momentum(maxOrder = momentum.Momentum.MAX_ORDER + 2)
        mmt = momentum.Momentum()
        self.assertTrue(math.isfinite(mmt._sMomentum[-1]))

    def testMaxOrder(self):
        mmt = momentum.Momentum()
        self.assertEqual(len(mmt._sMomentum) * 2, mmt.MAX_ORDER)



if __name__ == '__main__':
    unittest.main()