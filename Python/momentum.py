import datetime
import math
import os


class Normal:
    '''
    Calculate variance momentum.

    {BINDING_FACTOR} is default to 5 for leakage smaller than 5e-7
    {MAX_ORDER} is default to the max 122 before variance overflows

    The accurate calculation of higher-order momentum of is too slow, so approximate calculation is used.

    So far, manual run of unittest can only be carried out in VarianceArithemtic\Python folder.
    So both FILE_APPROX and FILE_PRECISE has two location
    '''
    BINDING_FACTOR:float = 5.0
    MAX_ORDER:int = 442
    FILE_APPROX = ('./Python/Output/NormalMomentum_5.txt', './Output/NormalMomentum_5.txt')
    FILE_PRECISE = ('./Python/NormalMomentum_5.txt', './NormalMomentum_5.txt')

    __slots__ = ('_sMomentum', '_maxOrder', '_binding', '_header', '_divid')

    def __init__(self, binding:float=BINDING_FACTOR, maxOrder:int=MAX_ORDER, divid:int=128,
                 readCached:bool=True) -> None:
        self._binding = binding
        self._maxOrder = maxOrder
        self._divid = divid
        self._header = f'n\tMomentum\t!!Diff\tSigma={self._binding}\n'

        self._sMomentum = []
        if readCached:
            if os.path.isfile(Normal.FILE_APPROX[0]):
                self._sMomentum = self._read(Normal.FILE_APPROX[0])
                if len(self._sMomentum) > maxOrder:
                    self._sMomentum = self._sMomentum[:maxOrder]
            elif os.path.isfile(Normal.FILE_APPROX[1]):
                self._sMomentum = self._read(Normal.FILE_APPROX[1])
                if len(self._sMomentum) > maxOrder:
                    self._sMomentum = self._sMomentum[:maxOrder]
        if not self._sMomentum:
            self._sMomentum = self.approxCalc()

        if readCached:
            if os.path.isfile(Normal.FILE_PRECISE[0]):
                sMomentum = self._read(Normal.FILE_PRECISE[0])
            elif os.path.isfile(Normal.FILE_PRECISE[1]):
                sMomentum = self._read(Normal.FILE_PRECISE[1])
            for i, mmt in enumerate(sMomentum):
                self._sMomentum[i] = mmt

    @property
    def binding(self):
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
                sMomentum = self._read(Normal.FILE_PRECISE[0])
            except:
                try:
                    sMomentum = self._read(Normal.FILE_PRECISE[1])
                except:
                    sMomentum = []
            x = sympy.symbols("x", is_real=True)
            while ((n := len(sMomentum)*2) < self.maxOrder):
                print(f'At {datetime.datetime.now()}, Calculate {n}/{self.maxOrder} momentum')
                density = 1/sympy.sqrt(2*sympy.pi) * sympy.exp(-x**2/2) * x**n
                mmt = sympy.integrate(density, (x, -self.binding, self.binding))
                if not math.isfinite(mmt):
                    print(f'The variance moment at {n} becomes inifinite')
                    break
                sMomentum.append(mmt.evalf())
                if self._write(Normal.FILE_PRECISE[0], sMomentum) or \
                   self._write(Normal.FILE_PRECISE[1], sMomentum):
                    return True
            print(f'Fail to write to either of {Normal.FILE_PRECISE}')
            return False
        except BaseException as ex:
            print(f'Exception when calculating momentum {sMomentum}: {ex}')
            return False

    def approxCalc(self) -> list[float]:       
        sMomentum = [0] * (self._maxOrder // 2)
        norm = 1.0 / math.sqrt(2*math.pi) / self._divid
        limit = int(math.ceil(self._divid * self.binding))
        for i in range(-limit, limit):
            x2 = (i + 0.5)**2 /self._divid**2
            pdf = norm * math.exp(- x2 * 0.5)
            sq = 1
            for j in range(self._maxOrder // 2):
                try:
                    sMomentum[j] += pdf * sq
                    sq *= x2
                    if not math.isfinite(sMomentum[j]):
                        raise ValueError(f'The {2*j}-th moment {sMomentum[j]} becomes overlow: {ex}')
                except BaseException as ex:
                    raise ValueError(f'The {2*j}-th moment {sMomentum[j]} becomes overlow: {ex}')
        if (not self._write(Normal.FILE_APPROX[0], sMomentum)) and \
           (not self._write(Normal.FILE_APPROX[1], sMomentum)):
            raise ValueError('Approx calc of momentum failed')
        return sMomentum


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



