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
        if os.path.isfile(momentum.Momentum.FILE_APPROX):
            os.remove(momentum.Momentum.FILE_APPROX)
        mmt = momentum.Momentum()
        self.assertEqual(len(mmt._sMomentum) * 2, mmt.MAX_ORDER)

    @unittest.skip('too slow')
    def testAnalyticCalc(self):
        '''
        n	Momentum	Uncertainty	!!Diff	Sigma=5.0
        0	0.999999426696856	-5.733031439580927e-07
        2	0.999984559501709	-1.5440498291052762e-05
        4	2.99958199862644	-0.00013933379118669098
        6	14.9886179961651	-0.000758800255658221
        8	104.688026048979	-0.00297118048591760
        10	936.384736336377	-0.00911668112552688
        12	10155.0446470894	-0.0230837280337300
        14	128385.894096893	-0.0499434336264284
        16	1835046.25357166	-0.0947096096142588
        18	28927232.3636749	-0.160542221361068
        20	492903566.233740	-0.247164078922660
        22	8933128674.00648	-0.350285337924554


        At 2024-08-09 23:10:11.278337, Calculate 6/126 momentum
        At 2024-08-09 23:10:16.347640, Calculate 8/126 momentum
        At 2024-08-09 23:10:40.965761, Calculate 10/126 momentum
        At 2024-08-09 23:12:03.707120, Calculate 12/126 momentum
        At 2024-08-09 23:19:46.447281, Calculate 14/126 momentum
        At 2024-08-09 23:42:27.935977, Calculate 16/126 momentum
        At 2024-08-10 00:12:49.837080, Calculate 18/126 momentum
        At 2024-08-10 02:50:42.099607, Calculate 20/126 momentum
        At 2024-08-10 06:57:59.741903, Calculate 22/126 momentum
        At 2024-08-10 18:41:27.228239, Calculate 24/122 momentum

        To be run by command line
        '''
        mmt = momentum.Momentum()
        self.assertTrue(mmt.analyticCalc())

if __name__ == '__main__':
    unittest.main()