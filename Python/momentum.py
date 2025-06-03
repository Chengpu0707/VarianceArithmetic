import datetime
import math
import os
import scipy.special
import scipy.stats

class Normal:
    '''
    Calculate variance momentum for the given {bounding}.
    Detect the {_maxOrder} for the bounding, which is 442 when {bounding}=5.

    The precise calculation of higher-order momentum of is too slow, 
    so approximate calculation is used with {_divid} divisions.
     *) {_filePrecise} stores the precise calculation.
     *) {_fileApprox} stores the approximate calculation.
     *) {_header} is the header of the output files.
    
    '''

    __slots__ = ('_sMomentum', '_maxOrder', '_binding', '_header', '_divid', '_fileApprox', '_filePrecise')

    @staticmethod
    def calcLow(bounding:float, n:int) -> float:
        if (n % 2) == 1:
            return 0.0
        if not (0 <= n <= 10):
            raise ValueError(f'Invalid n {n} for lower calculation')
        if not (0 <= bounding <= 10):
            raise ValueError(f'Invalid bounding {bounding} for lower calculation')
        first = scipy.special.erf(bounding / math.sqrt(2))
        dfrac = 1.0
        sum = 0
        num = bounding * 2 * scipy.stats.norm.pdf(bounding)
        for i in range(1, n // 2 + 1):
            dfrac *= 2*i - 1
            sum += num / dfrac
            num *= bounding**2
        return dfrac * (first - sum)
    

    def __init__(self, bounding:float=5.0, divid:int=128,
                 readCached:bool=True) -> None:
        self._binding = bounding
        self._divid = divid
        self._header = f'n\tMomentum\t!!Diff\tSigma={bounding}\n'
        self._fileApprox = f'./Python/Output/NormalMomentum_{bounding}.txt'
        self._filePrecise = f'./Python/NormalMomentum_{bounding}.txt'

        self._sMomentum = []
        if readCached:
            if os.path.isfile(self._fileApprox):
                self._sMomentum = self._read(self._fileApprox)
                self._maxOrder = len(self._sMomentum) * 2
        if not self._sMomentum:
            self._sMomentum = self.approxCalc()

        if readCached:
            if os.path.isfile(self._filePrecise):
                sMomentum = self._read(self._filePrecise)
                for i, mmt in enumerate(sMomentum):
                    self._sMomentum[i] = mmt

    @property
    def bounding(self):
        return self._binding
    
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


    def analyticCalc(self) -> bool:
        try:
            import sympy
            try:
                sMomentum = self._read(self._filePrecise)
            except:
                sMomentum = []
            x = sympy.symbols("x", is_real=True)
            while ((n := len(sMomentum)*2) < self.maxOrder):
                print(f'At {datetime.datetime.now()}, Calculate {n}/{self.maxOrder} momentum')
                density = 1/sympy.sqrt(2*sympy.pi) * sympy.exp(-x**2/2) * x**n
                mmt = sympy.integrate(density, (x, -self.bounding, self.bounding))
                if not math.isfinite(mmt):
                    print(f'The variance moment at {n} becomes inifinite')
                    break
                sMomentum.append(mmt.evalf())
                if self._write(self._filePrecise[0], sMomentum) or \
                   self._write(self._filePrecise[1], sMomentum):
                    return True
            print(f'Fail to write to either of {self._filePrecise}')
            return False
        except BaseException as ex:
            print(f'Exception when calculating momentum {sMomentum}: {ex}')
            return False

    def approxCalc(self) -> list[float]:
        norm = 1.0 / math.sqrt(2*math.pi) / self._divid
        limit = int(math.ceil(self._divid * self.bounding))
        self._maxOrder = 1000       
        sMomentum = [0] * (self._maxOrder // 2)
        for i in range(-limit, limit):
            x2 = (i + 0.5)**2 /self._divid**2
            pdf = norm * math.exp(- x2 * 0.5)
            sq = 1
            for j in range(self._maxOrder // 2):
                try:
                    sMomentum[j] += pdf * sq
                    sq *= x2
                    if not math.isfinite(sMomentum[j]):
                        self._maxOrder = j * 2
                except BaseException as ex:
                    self._maxOrder = j * 2
        if not self._write(self._fileApprox, sMomentum):
            raise ValueError('Approx calc of momentum failed')
        return sMomentum[:self.maxOrder // 2]


    def _write(self, file: str, sMomentum: list[float], withFac=False):
        try:
            with open(file, 'w') as f:
                f.write(self._header)
                dbFac = 1.
                for i, mmt in enumerate(sMomentum):
                    if not math.isfinite(mmt):
                        break
                    if withFac and math.isfinite(dbFac):
                        if i:
                            try:
                                dbFac *= 2*i - 1
                            except OverflowError:
                                pass
                        f.write(f'{i*2}\t{mmt}\t{mmt / dbFac}\n')
                    else:
                        f.write(f'{i*2}\t{mmt}\n')
                    f.flush()
                return True
        except BaseException as ex:
            return False

    def _read(self, file: str) -> list[float]:
        sMomentum = []
        with open(file) as f:
            header = next(f)
            if header != self._header:
                raise IOError(f'{file} has wrong header "{header}" != "{self._header}"')
            dbFac = 1.0
            for line in f:
                sWord = line.split('\t')
                i = int(sWord[0])
                if i != (len(sMomentum) * 2):
                    raise IOError(f'{file} has wrong index {i} vs {len(sMomentum)} at line "{line}"')
                momentum = float(sWord[1])
                match i:
                    case 0:
                        if not ((1 - 1e-06) < momentum < 1):
                            raise IOError(f'{file} has wrong momentum {momentum} <= {sMomentum[-1]} at index {i} of line "{line}"')
                    case 2:
                        if not ((1 - 2e-05) < momentum <= sMomentum[-1]):
                            raise IOError(f'{file} has wrong momentum {momentum} <= {sMomentum[-1]} at index {i} of line "{line}"')
                    case _:
                        if sMomentum[-1] >= momentum:
                            raise IOError(f'{file} has wrong momentum {momentum} <= {sMomentum[-1]} at index {i} of line "{line}"')
                        dbFac *= i - 1
                        if math.isfinite(dbFac):
                            if momentum >= dbFac:
                                raise IOError(f'{file} has wrong momentum {momentum} >= {dbFac} at index {i} of line "{line}"')
                sMomentum.append(momentum) 
        return sMomentum



