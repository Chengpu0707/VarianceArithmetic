"""Unit tests for recursiveSin.py — exercises RecursiveSin's recursive
generation of sin/cos values from sin(0)=0 and sin(pi/2)=1 and validates
self-consistency of summation identities.
"""
import unittest

from recursiveSin import RecursiveSin
from varDbl import VarDbl

class TestRecursiveSin (unittest.TestCase):

    def test_2(self):
        regr = RecursiveSin(2)
        regr.calc()

    def test_sum2(self):
        '''
        rSin = regr._sSin[1]
        rCos = regr._sSin[131071]
        Even though the uncertainty of rSin is much larger than the uncertainty of qSin,
            the result uncertainty for Sin ** 2 + Cos ** 2 is consistent
        '''
        rSin = VarDbl(0.0000119842258870683, 1.38098960702522E-12)
        rCos = VarDbl(0.999999999928189, 6.62008766811814E-17)
        qSin = VarDbl(0.0000119842249050697, 1.6940658945086E-21)
        qCos = VarDbl(0.999999999928189, 1.11022302462515E-16)
        rErr = rSin ** 2 + rCos ** 2 - 1
        qErr = qSin ** 2 + qCos ** 2 - 1
        self.assertAlmostEqual(rErr.value(), 0)
        self.assertAlmostEqual(qErr.value(), 0)
        self.assertAlmostEqual(rErr.uncertainty(), qErr.uncertainty())

    def test_consistency(self):
        regr4 = RecursiveSin(4)
        regr4.calc()
        regr5 = RecursiveSin(5)
        regr5.calc()
        for i in range(0, 5):
            self.assertAlmostEqual(regr4._sSin[i].value(), regr5._sSin[i << 1].value())
            self.assertAlmostEqual(regr4._sSin[i].uncertainty(), regr5._sSin[i << 1].uncertainty())






if __name__ == '__main__':
    unittest.main()



