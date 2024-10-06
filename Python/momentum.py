import datetime
import math
import os


class Momentum:
    '''
    Calculate variance momentum.

    {BINDING_FACTOR} is default to 5 for leakage smaller than 5e-7
    {MAX_ORDER} is default to the max 122 before variance overflows

    The accurate calculation of higher-order momentum of is too slow, so approximate calculation is used.

    So far, manual run of unittest can only be carried out in VarianceArithemtic\Python folder.
    So both FILE_APPROX and FILE_PRECISE has two location
    '''
    BINDING_FACTOR:float = 5.0
    MAX_ORDER:int = 244
    FILE_APPROX = ('./Python/Output/NormalMomentum_5.txt', './Output/NormalMomentum_5.txt')
    FILE_PRECISE = ('./Python/NormalMomentum_5.txt', './NormalMomentum_5.txt')

    __slots__ = ('_sMomentum', '_maxOrder', '_binding', '_header')

    def __init__(self, binding:float=BINDING_FACTOR, maxOrder:int=MAX_ORDER) -> None:
        self._binding = binding
        self._maxOrder = maxOrder
        self._header = f'n\tMomentum\tUncertainty\t!!Diff\tSigma={self._binding}\n'

        if os.path.isfile(Momentum.FILE_APPROX[0]):
            self._sMomentum = self._read(Momentum.FILE_APPROX[0])
        elif os.path.isfile(Momentum.FILE_APPROX[1]):
            self._sMomentum = self._read(Momentum.FILE_APPROX[1])
        else:
            self._sMomentum = self.approxCalc()
        if os.path.isfile(Momentum.FILE_PRECISE[0]):
            sMomentum = self._read(Momentum.FILE_PRECISE[0])
        elif os.path.isfile(Momentum.FILE_PRECISE[1]):
            sMomentum = self._read(Momentum.FILE_PRECISE[1])
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
                sMomentum = self._read(Momentum.FILE_PRECISE[0])
            except:
                try:
                    sMomentum = self._read(Momentum.FILE_PRECISE[1])
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
                if self._write(Momentum.FILE_PRECISE[0], sMomentum) or \
                   self._write(Momentum.FILE_PRECISE[1], sMomentum):
                    return True
            print(f'Fail to write to either of {Momentum.FILE_PRECISE}')
            return False
        except BaseException as ex:
            print(f'Exception when calculating momentum {sMomentum}: {ex}')
            return False

    def approxCalc(self) -> list[float]:       
        import varDbl
        sMomentum = []
        for j in range(self._maxOrder // 2):
            sMomentum.append(varDbl.VarDbl())
        divid = 64
        norm = varDbl.VarDbl(1.0 / math.sqrt(2*math.pi) / divid)
        limit = int(math.ceil(divid * self.binding))
        for i in range(-limit, limit + 1):
            x2 = i*i /divid /divid
            pdf = norm * varDbl.VarDbl(math.exp(- x2 * 0.5))
            sq = varDbl.VarDbl(1)
            for j in range(self._maxOrder // 2):
                try:
                    sMomentum[j] += pdf * sq
                    sq *= x2
                except BaseException as ex:
                    raise ValueError(f'The {2*j}-th moment {sMomentum[j]} becomes overlow: {ex}')
                if sMomentum[j].uncertainty() * self.binding > sMomentum[j].value():
                    raise ValueError(f'The {2*j}-th moment {sMomentum[j]} becomes unreliable')
        if (not self._write(Momentum.FILE_APPROX[0], sMomentum)) and \
           (not self._write(Momentum.FILE_APPROX[1], sMomentum)):
            raise ValueError('Approx calc of momentum failed')
        return [mmt.value() for mmt in sMomentum]


    def _write(self, file: str, sMomentum: list[float]):
        import varDbl
        try:
            with open(file, 'w') as f:
                f.write(self._header)
                dbFac = 1
                for i, mmt in enumerate(sMomentum):
                    unc = 0
                    if type(mmt) == varDbl.VarDbl:
                        unc = mmt.uncertainty()
                        mmt = mmt.value()
                    if i > 1:
                        dbFac *= 2*i - 1
                        if not (sMomentum[i - 1] < mmt < dbFac):
                            print(f'Invalid {i}-th momentum {mmt} vs !!={dbFac} and {i - 1}-th momentum {sMomentum[i - 1]}')
                            f.flush()
                            return False
                        if not math.isfinite(mmt):
                            break
                    f.write(f'{i*2}\t{mmt}\t{unc}\t{(mmt / dbFac) - 1}\n')
                f.flush()
                return True
        except:
            return False

    def _read(self, file: str) -> list[float]:
        sMomentum = []
        with open(file) as f:
            header = next(f)
            if header != self._header:
                raise IOError(f'{file} has wrong header "{header}" != "{self._header}"')
            dbFac = 1
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
                        dbFac *= i - 1
                        if not (sMomentum[-1] < momentum < dbFac):
                            raise IOError(f'{file} has wrong momentum {momentum} <= {sMomentum[-1]} at index {i} of line "{line}"')
                sMomentum.append(momentum) 
        return sMomentum



