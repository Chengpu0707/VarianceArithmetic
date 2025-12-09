import os
import unittest

import momentum
import varDbl

class TestNormal (unittest.TestCase):
    filePath = './Python/NormalMomentum_5.0.txt'

    @staticmethod
    def preciseNormalMomentum5():
        return ( 1*(1 - 5.733031e-07), 
                 1*(1 - 1.544050e-05),
                 3*(1 - 1.393338e-04),
                15*(1 - 7.588003e-04),
               105*(1 - 2.9711805e-03),
               945*(1 - 9.1166811e-03));   

    @unittest.skip('1 day slow')
    def testCalcPrecice(self):
        momentum.Normal.calcPreciseNorm(maxOrder=24)
        self.assertTrue(os.path.isfile(TestNormal.filePath))

    @unittest.skip('file may not exist')
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

    def testNormal_withoutVariance(self):
        mmt = momentum.Normal(withVariance=False)
        self.assertAlmostEqual(mmt.leakage, 5.7330314e-07)
        self.assertEqual(mmt.maxOrder, 448)
        sMomentum = TestNormal.preciseNormalMomentum5()
        for i, m in enumerate(sMomentum):
            self.assertAlmostEqual(mmt[i*2] / m, 1)
            self.assertEqual(mmt[i*2 + 1], 0)

    def testNormal_withVariance(self):
        mmt = momentum.Normal(withVariance=True)
        self.assertEqual(mmt.maxOrder, 250)
        sMomentum = TestNormal.preciseNormalMomentum5()
        for i, m in enumerate(sMomentum):
            self.assertAlmostEqual(mmt[i*2].value() / m, 1)
            self.assertEqual(mmt[i*2 + 1], 0)

    def testCompare(self):
        mmtV = momentum.Normal(withVariance=True)
        self.assertEqual(mmtV.maxOrder, 250)
        self.assertEqual(mmtV.bounding, 5)
        mmtF = momentum.Normal(withVariance=False)
        self.assertGreater(mmtF.maxOrder, 250)
        self.assertEqual(mmtF.bounding, 5)
        sMomentum =TestNormal.preciseNormalMomentum5()
        with open('./Python/Output/NormalMomentum_compare.txt', 'w') as f:
            f.write('Order\tVar Value\tVar Uncertainty\tVar Precison\tVar Diff\tFloat Diff\tFloat Uncertainty\tFloat Precison\n')
            for n in range(0, mmtV.maxOrder, 2):
                unc = varDbl.VarDbl(mmtF[n]).uncertainty()
                f.write(f'{n}\t{mmtV[n].value()}\t{mmtV[n].uncertainty()}\t{mmtV[n].precision()}'
                        f'\t{mmtV[n].value() - sMomentum[n >> 1] if n < len(sMomentum) * 2 else ""}'
                        f'\t{mmtF[n] - mmtV[n].value()}\t{unc}\t{unc/mmtF[n]}\n')


class TestUniform (unittest.TestCase):

    def testUniform(self):
        mmt = momentum.Uniform()
        self.assertEqual(mmt.maxOrder, 653)
        self.assertEqual(mmt[0], 1)
        self.assertEqual(mmt[1], 0)
        self.assertEqual(mmt[2], 3**1 /3)
        self.assertEqual(mmt[3], 0)
        self.assertEqual(mmt[4], 3**2 /5)
        self.assertEqual(mmt[5], 0)
        self.assertEqual(mmt[6], 3**3 /7)
        self.assertEqual(mmt[7], 0)
        self.assertEqual(mmt[8], 3**4 /9)


if __name__ == '__main__':
    unittest.main()