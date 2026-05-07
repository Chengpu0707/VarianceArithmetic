"""Bounded moments for symmetric distributions used by Taylor expansion:
abstract Moment base class with NormalMoment and Uniform implementations
returning the moment table truncated by the bounding factor.
"""
import abc
import datetime
import math
import os
import scipy.special
import scipy.stats
import sympy

from indexSin import OUTDIR
import varDbl

class Moment (abc.ABC):

    @property
    @abc.abstractmethod
    def bounding(self):
        pass

    @property
    @abc.abstractmethod
    def leakage(self):
        pass

    @property
    @abc.abstractmethod
    def maxOrder(self):
        pass

    @abc.abstractmethod
    def __getitem__(self, n:int) -> float:
        pass


class Normal (Moment):
    '''
    Calculate variance moment for the given {bounding}.
    Detect the {_maxOrder} for the bounding, which is 448 when {bounding}=5.
    '''
    __slots__ = ('_sMoment', '_maxOrder', '_bounding')

    @staticmethod
    def readPreciseNorm(filePath:str):
        sMoment = []
        with open(filePath) as f:
            hdr = next(f)
            sHdr = hdr.split('\t')
            if (len(sHdr) != 4) or (sHdr[0] != 'n') or (sHdr[1] != 'Moment') or (sHdr[2] != 'Bounding:'):
                raise ValueError(f'Invalid header "{hdr}" in {filePath}')
            bounding = float(sHdr[3])
            if not (0 < bounding <= 8):
                raise ValueError(f'Invalid bounding {bounding} in header "{hdr}" in {filePath}')
            for i, line in enumerate(f):
                n, mmt = map(float, line.split('\t'))
                if n != i * 2:
                    raise ValueError(f'Invalid index {n} at line {i+2} in {filePath}, expected {i*2}: {line}')
                sMoment.append(mmt)
        return bounding, sMoment

    @staticmethod
    def calcPreciseNorm(bounding:float=5.0, maxOrder:int=1000000):
        filePath = f'{OUTDIR}/Python/NormalMoment_{bounding}.txt'
        HEADER = f'n\tMoment\tBounding={bounding}'
        b, sMoment = Normal.readPreciseNorm(filePath)
        if b != bounding:
            sMoment = []
        x = sympy.symbols("x", is_real=True)
        # Normalized per Formula (2.2): ζ(n, κ) = ∫z^n ρ dz / ∫ρ dz.
        density_const = 1/sympy.sqrt(2*sympy.pi) * sympy.exp(-x**2/2)
        norm = sympy.integrate(density_const, (x, -bounding, bounding)).evalf()
        while ((n := len(sMoment)*2) < maxOrder):
            print(f'At {datetime.datetime.now()}, Calculate {n}/{maxOrder} moment')
            mmt = sympy.integrate(density_const * x**n, (x, - bounding, bounding))
            if not math.isfinite(mmt):
                break
            sMoment.append((mmt / norm).evalf())
            with open(filePath, 'w') as f:
                f.write(HEADER)
                for i, mmt in enumerate(sMoment):
                    f.write(f'{i*2}\t{mmt}\n')

    def __init__(self, bounding:float=5, maxOrder:int=1000000, withVariance:bool=False):
        self._bounding = bounding
        filePath = f'{OUTDIR}/Python/Output/NormalMoment_{bounding}_{"var" if withVariance else "float"}.txt'
        # Header bumped from 'Value' to 'NormalizedValue' to invalidate caches
        # written before Formula (2.2) was normalized.
        HEADER = 'Order\tNormalizedValue\tUncertainty\n'
        if os.path.isfile(filePath):
            try:
                with open(filePath) as f:
                    hdr = next(f)
                    if hdr != HEADER:
                        raise NotImplementedError(f'Invalid header {hdr} vs {HEADER}')
                    n = -1
                    prevVal = 0
                    prevUnc = 0
                    sMoment = []
                    for line in f.readlines():
                        n += 1
                        nn, val, unc = map(float, line.strip().split('\t'))
                        if nn != n:
                            raise NotImplementedError(f'Invalid index {nn} vs {n}')
                        if (n & 1) == 0:
                            if val <= 0:
                                raise NotImplementedError(f'Invalid value {val} for index={n}')
                            if (unc <= 0) or (val/bounding < unc):
                                raise NotImplementedError(f'Invalid {val}+/-{unc} for index={n}')
                            if (val <= prevVal) or (unc <= prevUnc):
                                raise NotImplementedError(f'Invalid {val}+/-{unc} vs {prevVal}+/-{prevUnc} for index={n}')
                            if withVariance:
                                sMoment.append(varDbl.VarDbl(val, unc))
                            else:
                                sMoment.append(val)
                        else:
                            if val != 0:
                                raise NotImplementedError(f'Invalid value {val} for index={n}')
                    self._maxOrder = len(sMoment) * 2
                    self._sMoment = sMoment
                    return
            except BaseException as ex:
                os.remove(filePath)

        term = 2 * scipy.stats.norm.pdf(bounding) * self._bounding
        bounding2 = self._bounding**2
        if withVariance:
            term = varDbl.VarDbl(term)
            bounding2 = varDbl.VarDbl(bounding2)
        sTerm = []
        for n in range(maxOrder):
            try:
                sTerm.append(term * (1/(2*n + 1)))
                if withVariance:
                    if (not math.isfinite(sTerm[-1].value())) or (not math.isfinite(sTerm[-1].variance())):
                        del sTerm[-1]
                        break
                else:
                    if not math.isfinite(sTerm[-1]):
                        del sTerm[-1]
                        break
                term *= bounding2
            except:
                break
        self._maxOrder = n * 2
        self._sMoment = sTerm[:]
        for j in range(2, maxOrder):
            for i in range(n):
                sTerm[i] *= 1/(2*i - 1 + 2*j) * bounding2
                prev = self._sMoment[i]
                self._sMoment[i] += sTerm[i]
                if (prev.value() == self._sMoment[i].value()) if withVariance else (prev == self._sMoment[i]):
                    n = i
                    break
            if n == 0:
                break
        # Normalize per Formula (2.2): divide unnormalized moments by ∫ρ dz = 1 - leakage.
        norm_factor = 1 - self.leakage
        for i in range(len(self._sMoment)):
            self._sMoment[i] = self._sMoment[i] * (1 / norm_factor)
        with open(filePath, 'w') as f:
            f.write(HEADER)
            for n in range(self.maxOrder):
                mmt = self[n]
                if type(mmt) == varDbl.VarDbl:
                    f.write(f'{n}\t{mmt.value()}\t{mmt.uncertainty()}\n')
                else:
                    f.write(f'{n}\t{mmt}\t{math.ulp(mmt)}\n')

    @property
    def bounding(self):
        return self._bounding

    @property
    def leakage(self):
        return 1 - scipy.special.erf(self.bounding/math.sqrt(2))

    @property
    def maxOrder(self):
        return self._maxOrder

    def __getitem__(self, n:int) -> float:
        if n < 0:
            return IndexError()
        if (n % 2) == 1:
            return 0
        n //= 2
        if n >= len(self._sMoment):
            return IndexError()
        return self._sMoment[n]


class Uniform:
    '''
    Pre-calculated variance moment for uniform distribution [-1, 1].
    '''
    __slots__ = ('_bounding', '_sMoment', '_maxOrder')

    def __init__(self, bounding=math.sqrt(3)):
        self._bounding = bounding
        self._sMoment = []
        fac = 1
        for n in range(10000):
            try:
                mmt = fac/(2*n + 1)
            except OverflowError:
                break
            if not math.isfinite(mmt):
                break
            self._sMoment.append(mmt)
            fac *= bounding**2
        self._maxOrder = len(self._sMoment)

    @property
    def bounding(self):
        return self._bounding

    @property
    def leakage(self):
        return 1 - self.bounding

    @property
    def maxOrder(self):
        return self._maxOrder

    def __getitem__(self, n:int) -> float:
        if n < 0 or n >= self._maxOrder or (n % 2) == 1:
            return 0
        return self._sMoment[n >> 1]


NORMAL = Normal(bounding=5.0)

UNIFORM = Uniform()
