import unittest

from momentum import Momentum
from varDbl import UncertaintyException, validate

class TestMomuntum (unittest.TestCase):
    momentum = Momentum()

    def test_max_order(self):
        try:
            Momentum(maxOrder=Momentum.MAX_ORDER+1)
        except UncertaintyException:
            pass
        except BaseException as ex:
            raise ex

    def test_even(self):
        validate(self, TestMomuntum.momentum.factor(0),    1 * 0.99999947194737793)
        validate(self, TestMomuntum.momentum.factor(2),    1 * 0.99998569318184616)
        validate(self, TestMomuntum.momentum.factor(4),    3 * 0.99987013368935596)
        validate(self, TestMomuntum.momentum.factor(6),   15 * 0.99928863789470035)
        validate(self, TestMomuntum.momentum.factor(8),  105 * 0.99719860134891214)
        validate(self, TestMomuntum.momentum.factor(10), 945 * 0.99135593485973217) 

    def test_odd(self):
        validate(self, TestMomuntum.momentum.factor(1), 0)
        validate(self, TestMomuntum.momentum.factor(3), 0)
        validate(self, TestMomuntum.momentum.factor(5), 0)
        validate(self, TestMomuntum.momentum.factor(7), 0)
        validate(self, TestMomuntum.momentum.factor(9), 0)

if __name__ == '__main__':
    unittest.main()