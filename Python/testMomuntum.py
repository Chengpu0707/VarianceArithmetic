import unittest

import momentum

class TestMomuntum (unittest.TestCase):

    def testEven(self):
        _momentum = momentum.Momentum(200, 5)
        self.assertEquals(  1 * 0.99999947194737793, _momentum.factor(0))
        self.assertEquals(  1 * 0.99998569318184616, _momentum.factor(2))
        self.assertEquals(  3 * 0.99987013368935596, _momentum.factor(4))
        self.assertEquals( 15 * 0.99928863789470035, _momentum.factor(6))
        self.assertEquals(105 * 0.99719860134891214, _momentum.factor(8))
        self.assertEquals(945 * 0.99135593485973217, _momentum.factor(10)) 

    def testOdd(self):
        _momentum = momentum.Momentum(200, 5)
        self.assertEquals(0, _momentum.factor(1))
        self.assertEquals(0, _momentum.factor(3))
        self.assertEquals(0, _momentum.factor(5))
        self.assertEquals(0, _momentum.factor(7))
        self.assertEquals(0, _momentum.factor(9))

if __name__ == '__main__':
    unittest.main()