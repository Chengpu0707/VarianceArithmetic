"""Unit tests for moment.py — verifies NormalMoment values against
high-precision reference table, checks bounding/leakage, and writes the
generated moment table for cross-implementation comparison.
"""
import os
import unittest

from indexSin import OUTDIR
import moment
import varDbl


class TestNormal (unittest.TestCase):
    filePath = f'{OUTDIR}/Python/NormalMoment_5.0.txt'

    @staticmethod
    def preciseNormalMoment5():
        # Normalized per Formula (2.2): deficits relative to the ideal
        # Gaussian even moments; ζ(0,κ)=1 by construction.
        return ( 1*(1 - 0),
                 1*(1 - 1.4867204e-05),
                 3*(1 - 1.3876057e-04),
                15*(1 - 7.5822739e-04),
               105*(1 - 2.9706089e-03),
               945*(1 - 9.1161130e-03))

    @unittest.skip('1 day slow')
    def testCalcPrecice(self):
        moment.Normal.calcPreciseNorm(maxOrder=24)
        self.assertTrue(os.path.isfile(TestNormal.filePath))

    @unittest.skip('file may not exist')
    def testReadPrecice(self):
        bounding, sMoment = moment.Normal.readPreciseNorm(TestNormal.filePath)
        self.assertEqual(bounding, 5)
        self.assertEqual(len(sMoment), 13)
        self.assertAlmostEqual(sMoment[0],   1*(1 - 0))
        self.assertAlmostEqual(sMoment[1],   1*(1 - 1.4867204e-05))
        self.assertAlmostEqual(sMoment[2],   3*(1 - 1.3876057e-04))
        self.assertAlmostEqual(sMoment[3],  15*(1 - 7.5822739e-04))
        self.assertAlmostEqual(sMoment[4], 105*(1 - 2.9706089e-03))
        self.assertAlmostEqual(sMoment[5], 945*(1 - 9.1161130e-03))

    def testNormal_withoutVariance(self):
        mmt = moment.Normal(withVariance=False)
        self.assertAlmostEqual(mmt.leakage, 5.7330314e-07)
        self.assertEqual(mmt.maxOrder, 448)
        sMoment = TestNormal.preciseNormalMoment5()
        for i, m in enumerate(sMoment):
            self.assertAlmostEqual(mmt[i*2] / m, 1)
            self.assertEqual(mmt[i*2 + 1], 0)

    def testNormal_withVariance(self):
        mmt = moment.Normal(withVariance=True)
        self.assertEqual(mmt.maxOrder, 250)
        sMoment = TestNormal.preciseNormalMoment5()
        for i, m in enumerate(sMoment):
            self.assertAlmostEqual(mmt[i*2].value() / m, 1)
            self.assertEqual(mmt[i*2 + 1], 0)

    def testCompare(self):
        mmtV = moment.Normal(withVariance=True)
        self.assertEqual(mmtV.maxOrder, 250)
        self.assertEqual(mmtV.bounding, 5)
        mmtF = moment.Normal(withVariance=False)
        self.assertGreater(mmtF.maxOrder, 250)
        self.assertEqual(mmtF.bounding, 5)
        sMoment = TestNormal.preciseNormalMoment5()
        with open(f'{OUTDIR}/Python/Output/NormalMoment_compare.txt', 'w') as f:
            f.write('Order\tVar Value\tVar Uncertainty\tVar Precison\tVar Diff\tFloat Diff\tFloat Uncertainty\tFloat Precison\n')
            for n in range(0, mmtV.maxOrder, 2):
                unc = varDbl.VarDbl(mmtF[n]).uncertainty()
                f.write(f'{n}\t{mmtV[n].value()}\t{mmtV[n].uncertainty()}\t{mmtV[n].precision()}'
                        f'\t{mmtV[n].value() - sMoment[n >> 1] if n < len(sMoment) * 2 else ""}'
                        f'\t{mmtF[n] - mmtV[n].value()}\t{unc}\t{unc/mmtF[n]}\n')


class TestUniform (unittest.TestCase):

    def testUniform(self):
        mmt = moment.Uniform()
        self.assertEqual(mmt.maxOrder, 647)
        self.assertAlmostEqual(mmt[0], 1, places=13)
        self.assertEqual(mmt[1], 0)
        self.assertAlmostEqual(mmt[2], 3**1 /3, places=13)
        self.assertEqual(mmt[3], 0)
        self.assertAlmostEqual(mmt[4], 3**2 /5, places=13)
        self.assertEqual(mmt[5], 0)
        self.assertAlmostEqual(mmt[6], 3**3 /7, places=13)
        self.assertEqual(mmt[7], 0)
        self.assertAlmostEqual(mmt[8], 3**4 /9, places=13)


if __name__ == '__main__':
    unittest.main()
