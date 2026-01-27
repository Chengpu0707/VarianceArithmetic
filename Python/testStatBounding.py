'''
Sample from normal distribution, and calculate bounding leakage.
'''
import functools
import math
import os
import unittest

from indexSin import OUTDIR
import momentum
import taylor
import varDbl


class TestBoundingFactor (unittest.TestCase):

    def assert_func(self, fileName:str):
        sMomentum = {}
        filePath = f'{OUTDIR}/Python/Output/{fileName}Leakage.txt'
        match fileName:
            case 'Normal':
                with open(filePath) as f:
                    hdr = next(f)
                    self.assertEqual(hdr, 'Samples\tKappa\tCount\tMean\tDeviation\tRange\n')
                    for line in f:
                        n, kappa, cnt, leak, std, bounding = map(float, line.split('\t'))
                        if kappa != 5.0:
                            continue
                        samples = int(n)
                        if samples in sMomentum:
                            raise ValueError(f'Duplicate samples {samples}')
                        sMomentum[samples] = momentum.Normal(bounding=bounding)
            case 'Uniform':
                with open(f'{OUTDIR}/Python/Output/UniformLeakage.txt') as f:
                    hdr = next(f)
                    self.assertEqual(hdr, 'Samples\tCount\tMean\tDeviation\n')
                    for line in f:
                        n, cnt, leak, std = map(float, line.split('\t'))
                        samples = int(n)
                        if samples in sMomentum:
                            raise ValueError(f'Duplicate samples {samples}')
                        bounding = 1 - leak
                        sMomentum[samples] = momentum.Uniform(bounding=bounding)
            case _:
                raise ValueError(f'Invalid fileName {fileName}')

        sFunc = {
            'x': lambda mmt, var: taylor.Taylor.polynominal1d(var, (0,1), momentum=mmt),
            'sin(x)': lambda mmt, var: taylor.Taylor.sin(var, momentum=mmt),
            'exp(x)': lambda mmt, var: taylor.Taylor.exp(var, momentum=mmt), 
            'log(x)': lambda mmt, var: taylor.Taylor.log(var, momentum=mmt), 
        }
        var = varDbl.VarDbl(1, 0.1)
        zero = varDbl.VarDbl(1, 0)

        def powFunc(exp, mmt, var):
            return taylor.Taylor.pow(var, exp, momentum=mmt)

        for exp in (2, 0.5, -1, -2):
            sFunc[f'x^{exp}'] = functools.partial(powFunc, exp)

        filePath = f'{OUTDIR}/Python/Output/{fileName}Bounding.txt'
        with open(filePath, 'w') as f:
            f.write('Samples\tBounding\tLeakage\tFunction\tStable Bias\tStable Variance'
                    '\tOutput Bias\tOutput Variance\tBias Ratio\tVariance Ratio'
                    '\tAdjust\tAdjusted Bias\tAdjusted Variance\tNormalized Bias\tNormalized Variance\n')
            for calc, func in sFunc.items():
                val = func(momentum.NORMAL, zero)
                match fileName:
                    case 'Normal':
                        stable = func(momentum.NORMAL, var)
                    case 'Uniform':
                        stable = func(momentum.UNIFORM, var)
                bias = stable.value() - val.value()
                for samples, mmt in sMomentum.items():
                    res = func(mmt, var)
                    adj = (res.value() - val.value()) / mmt[2]
                    adjVar = res.variance() / mmt[2]
                    if bias == 0:
                        f.write(f'{samples}\t{mmt.bounding}\t{mmt.leakage}\t{calc}\t{bias}\t{stable.variance()}'
                                f'\t{res.value() - val.value()}\t{res.variance()}\t1\t{res.variance()/stable.variance()}'
                                f'\t{mmt[2]}\t{adj}\t{adjVar}\t1\t{adjVar/stable.variance()}\n')
                    else:
                        f.write(f'{samples}\t{mmt.bounding}\t{mmt.leakage}\t{calc}\t{bias}\t{stable.variance()}'
                                f'\t{res.value() - val.value()}\t{res.variance()}\t{(res.value() - val.value())/bias}\t{res.variance()/stable.variance()}'
                                f'\t{mmt[2]}\t{adj}\t{adjVar}\t{adj/bias}\t{adjVar/stable.variance()}\n')
        self.assertTrue(os.path.isfile(filePath))

    def test_Normal(self):
        self.assert_func('Normal')

    def test_Uniform(self):
        self.assert_func('Uniform')

    def test_Uniform_Adj(self):
        mmt = momentum.Uniform(0.5)
        ideal = taylor.Taylor.pow(varDbl.VarDbl(1, 0.1), -2, momentum=momentum.UNIFORM, dumpPath=f'{OUTDIR}/Python/Output/pow_1_0.1_-2.txt')
        half = taylor.Taylor.pow(varDbl.VarDbl(1, 0.1), -2, momentum=mmt, dumpPath=f'{OUTDIR}/Python/Output/pow_1_0.1_-2_mmt_0.5.txt')
        self.assertAlmostEqual((half.value() - 1) / mmt[2] / (ideal.value() - 1), 1, delta=0.03)
        self.assertAlmostEqual(half.variance() / mmt[2] / ideal.variance(), 1, delta=0.07)

        


if __name__ == '__main__':
    unittest.main()
