import datetime
import math
import numpy
import scipy.special
import scipy.stats
import sympy

class Normal:
    '''
    Calculate variance momentum for the given {bounding}.
    Detect the {_maxOrder} for the bounding, which is 448 when {bounding}=5.   
    '''

    __slots__ = ('_sMomentum', '_maxOrder', '_bounding')

    @staticmethod
    def getPath(bounding:float):
        return f'./Python/NormalMomentum_{bounding}.txt'

    @staticmethod
    def readPreciseNorm(filePath:str):
        sMomentum = []
        with open(filePath) as f:
            hdr = next(f)
            sHdr = hdr.split('\t')
            if (len(sHdr) != 4) or (sHdr[0] != 'n') or (sHdr[1] != 'Momentum') or (sHdr[2] != 'Bounding:'):
                raise ValueError(f'Invalid header "{hdr}" in {filePath}')
            bounding = float(sHdr[3])
            if not (0 < bounding <= 8):
                raise ValueError(f'Invalid bounding {bounding} in header "{hdr}" in {filePath}')
            for i, line in enumerate(f):
                n, mmt = map(float, line.split('\t'))
                if n != i * 2:
                    raise ValueError(f'Invalid index {n} at line {i+2} in {filePath}, expected {i*2}: {line}')
                sMomentum.append(mmt)
        return bounding, sMomentum

    @staticmethod
    def calcPreciseNorm(bounding:float=5.0, maxOrder:int=1000000):
        filePath = Normal.getPath(bounding)
        HEADER = f'n\tMomentum\tBounding={bounding}'
        b, sMomentum = Normal.readPreciseNorm(filePath)
        if b != bounding:
            sMomentum = []
        x = sympy.symbols("x", is_real=True)
        while ((n := len(sMomentum)*2) < maxOrder):
            print(f'At {datetime.datetime.now()}, Calculate {n}/{maxOrder} momentum')
            density = 1/sympy.sqrt(2*sympy.pi) * sympy.exp(-x**2/2) * x**n
            mmt = sympy.integrate(density, (x, - bounding, bounding))
            if not math.isfinite(mmt):
                break
            sMomentum.append(mmt.evalf())
            with open(filePath, 'w') as f:
                f.write(HEADER)
                for i, mmt in enumerate(sMomentum):
                    f.write(f'{i*2}\t{mmt}\n')

    @staticmethod
    def calcLow(bounding:float, n:int, readCached:bool=True) -> float:
        if (n % 2) == 1:
            return 0.0
        if not (0 <= n <= 20):
            raise ValueError(f'Invalid n {n} for lower calculation')
        if not (0 <= bounding <= 8):
            raise ValueError(f'Invalid bounding {bounding} for lower calculation')
        if readCached:
            try:
                b, sMomentum = Normal.readPreciseNorm(Normal.getPath(bounding))
                if (b == bounding) and (n < len(sMomentum)*2):
                    return sMomentum[n // 2]
            except:
                pass
        first = scipy.special.erf(bounding / math.sqrt(2))
        dfrac = 1.0
        sum = 0
        num =  2 * scipy.stats.norm.pdf(bounding) * bounding
        for i in range(1, n // 2 + 1):
            dfrac *= 2*i - 1
            sum += num / dfrac
            num *= bounding**2
        return dfrac * (first - sum)
    

    def __init__(self, bounding:float=5.0, maxOrder:int=1000000) -> None:
        self._bounding = numpy.float64(bounding)
        term = 2 * scipy.stats.norm.pdf(bounding) * self._bounding
        bounding2 = self._bounding**2
        sTerm = []
        for n in range(maxOrder):
            try:
                sTerm.append(term /(2*n + 1))
                if not math.isfinite(sTerm[-1]):
                    del sTerm[-1]
                    break
                term *= bounding2
            except:
                break
        self._maxOrder = n * 2
        self._sMomentum = sTerm[:]
        for j in range(2, maxOrder):
            for i in range(n):
                sTerm[i] = sTerm[i] / (2*i - 1 + 2*j) * bounding2
                prev = self._sMomentum[i]
                self._sMomentum[i] += sTerm[i]
                if prev == self._sMomentum[i]:
                    n = i
                    break
            if n == 0:
                break


    @property
    def bounding(self):
        return self._bounding
    
    @property
    def maxOrder(self):
        return self._maxOrder

    def __getitem__(self, n:int) -> float:
        if n < 0:
            return IndexError()
        if (n % 2) == 1:
            return 0
        n //= 2
        if n >= len(self._sMomentum):
            return IndexError()
        return self._sMomentum[n]


