import unittest

from momentum import Momentum
from varDbl import VarianceException, validate

class TestMomuntum (unittest.TestCase):
    momentum = Momentum()

    def test_max_order(self):
        try:
            Momentum(maxOrder=Momentum.MAX_ORDER+1)
        except VarianceException:
            pass
        except BaseException as ex:
            raise ex

    def test_even(self):
        validate(self, TestMomuntum.momentum.factor(0),    1 * 0.99999947194737793, 1.0631430356312753e-17)
        validate(self, TestMomuntum.momentum.factor(2),    1 * 0.99998569318184616, 1.0096418147444105e-17)
        validate(self, TestMomuntum.momentum.factor(4),    3 * 0.99987013368935596, 2.998915588392244e-17)
        validate(self, TestMomuntum.momentum.factor(6),   15 * 0.99928863789470035, 1.4926804198403977e-16)
        validate(self, TestMomuntum.momentum.factor(8),  105 * 0.99719860134891214, 1.0409495599106208e-15)
        validate(self, TestMomuntum.momentum.factor(10), 945 * 0.99135593485973217, 9.33728878727979e-15) 

    def test_odd(self):
        validate(self, TestMomuntum.momentum.factor(1), 0)
        validate(self, TestMomuntum.momentum.factor(3), 0)
        validate(self, TestMomuntum.momentum.factor(5), 0)
        validate(self, TestMomuntum.momentum.factor(7), 0)
        validate(self, TestMomuntum.momentum.factor(9), 0)

    def test_dump(self):
        with open('./Python/Output/Momentum.txt', 'w') as f:
            f.write('2n\t(2n-1)!!\tValue\tUncertainty\n')
            var = TestMomuntum.momentum.factor(0)
            f.write(f'0\t1\t{var.value()}\t{var.uncertainty()}\n')
            fac = 1.0
            for i in range(2, Momentum.MAX_ORDER, 2):
                fac *= (i - 1)
                var = TestMomuntum.momentum.factor(i)
                f.write(f'{i}\t{fac}\t{var.value()}\t{var.uncertainty()}\n')

if __name__ == '__main__':
    unittest.main()