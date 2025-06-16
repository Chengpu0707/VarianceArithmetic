import math
import os
import unittest

import momentum

class TestNormal (unittest.TestCase):
    filePath = './Python/NormalMomentum_5.0.txt'

    @unittest.skip('1 day slow')
    def testCalcPrecice(self):
        momentum.Normal.calcPreciseNorm(maxOrder=24)
        self.assertTrue(os.path.isfile(TestNormal.filePath))

    def testReadPrecice(self):
        bounding, sMomentum = momentum.Normal.readPreciseNorm(TestNormal.filePath)
        self.assertEqual(bounding, 5)
        self.assertEqual(len(sMomentum), 13)
        self.assertAlmostEqual(sMomentum[0],    1*(1 - 5.733031e-07))
        self.assertAlmostEqual(sMomentum[1],    1*(1 - 1.544050e-05))
        self.assertAlmostEqual(sMomentum[2],    3*(1 - 1.393338e-04))
        self.assertAlmostEqual(sMomentum[3],   15*(1 - 7.588003e-04))
        self.assertAlmostEqual(sMomentum[4],  105*(1 - 2.9711805e-03))
        self.assertAlmostEqual(sMomentum[5],  945*(1 - 9.1166811e-03))

    def testNormal(self):
        mmt = momentum.Normal()
        self.assertEqual(mmt.maxOrder, 448)
        bounding, sMomentum = momentum.Normal.readPreciseNorm(TestNormal.filePath)
        self.assertEqual(bounding, 5)
        for i, m in enumerate(sMomentum):
            self.assertAlmostEqual(mmt[i*2] / m, 1)
            self.assertEqual(mmt[i*2 + 1], 0)


    def testCalcLow(self):
        bounding, sMomentum = momentum.Normal.readPreciseNorm(TestNormal.filePath)
        self.assertEqual(bounding, 5)
        for i in range(11):
            self.assertAlmostEqual(momentum.Normal.calcLow(5, i*2, False) / sMomentum[i], 1)
            self.assertEqual(momentum.Normal.calcLow(5, i*2 + 1), 0)



if __name__ == '__main__':
    unittest.main()