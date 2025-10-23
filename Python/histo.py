
import math
import typing

class Stat:
    __slots__ = ['_count', '_min', '_max', '_minAt', '_maxAt', '_sum', '_sum2']

    def __init__(self) -> None:
        self._count = 0
        self._min = 0.0
        self._minAt = None
        self._max = 0.0
        self._maxAt = None
        self._sum = 0.0
        self._sum2 = 0.0

    def __str__(self) -> str:
        match self._count:
            case 0:
                return '"Stat: 0"'
            case 1:
                return f'"Stat: 1, {self._sum}"'
            case _:
                return f'"Stat: {self._count}, {self.mean()}+/-{self.dev()}"'
            
    def __add__(self, other):
        if not isinstance(other, Stat):
            return NotImplemented(f'{self} + {type(other)}: {other}')
        ret = Stat()
        ret._count = self._count + other._count
        if self._min > other._min:
            self._min = other._min
            self._minAt = other._minAt
        elif self._min == other._min:
            if (self._minAt is None) and (other._minAt is not None):
                self._minAt = other._minAt
        if self.max < other.max:
            self._max = other._max
            self._maxAt = other._maxAt
        elif self._max == other._max:
            if (self._maxAt is None) and (other._maxAt is not None):
                self._maxAt = other._maxAt
        ret._sum = self._sum + other._sum
        ret._sum2 = self._sum2 + other._sum2
        return ret

    def accum(self, value:float, at=None) ->bool:
        value = float(value)
        if not math.isfinite(value):
            return False
        if (not self._count) or (self._min > value):
            self._min = value
            if at is not None:
                self._minAt = at    
        if (not self._count) or (self._max < value):
            self._max = value
            if at is not None:
                self._maxAt = at
        self._count += 1
        self._sum += value
        self._sum2 += value * value
        return True

    def count(self):
        return self._count
    
    def min(self):
        return self._min if self._count else float('nan')
    
    def max(self):
        return self._max if self._count else float('nan')
    
    def minAt(self):
        return self._minAt
    
    def maxAt(self):
        return self._maxAt

    def mean(self):
        return self._sum / self._count if self._count else float('nan')

    def dev(self):
        if not self._count:
            return float('nan')
        var = self._sum2 / self._count - self.mean()**2
        return math.sqrt(var) if var >= 0 else 0.0        

class Histo:
    '''
    A histogram container specific for normalized error.
    The symmetric range is "devs".
    Within each dev, the bucket count is "divids"
    '''
    __slots__ = ['_stat', '_divids', '_devs', '_half', '_less', '_more', '_sHisto']

    def __init__(self, divids:int, devs:float) -> None:
        self._stat = Stat()
        self._divids = abs(int(divids))
        self._half = round(self._divids * abs(devs))
        self._devs = self._half / self._divids
        self._less = 0
        self._more = 0
        self._sHisto = [0] *(2*self._half + 1)


    def __str__(self) -> str:
        return str(self._stat)
    
    def __add__(self, other):
        if not isinstance(other, Histo):
            return NotImplemented(f'{self} + {type(other)}: {other}')
        if self._divids != other._divids:
            raise ValueError(f'Histo.divids {self._divids} + {other._divids}')
        if self._devs != other._devs:
            raise ValueError(f'Histo.devs {self._devs} + {other._devs}')
        if len(self._sHisto) != len(other._sHisto):
            raise ValueError(f'Histo.len {len(self._sHisto)} + {other._sHisto}')
        ret = Histo(self._divids, self._half /self._divids)
        ret.stat = self.stat + other.stat
        ret.less = self.less + other.less
        ret.more = self.more + other.more
        ret._sHisto = [it1 + it2 for it1, it2 in zip(self._sHisto, other._sHisto)]
        return ret
    
    def range(self):
        return self._devs
    
    def less(self):
        return self._less

    def more(self):
        return self._more

    def stat(self):
        return self._stat
    
    def accum(self, value:float, at=None) ->bool:
        if not self._stat.accum(value, at):
            return False
        if value < -self._devs:
            self._less += 1
        elif value > self._devs:
            self._more += 1
        else:
            self._sHisto[self._half + round(value*self._divids)] += 1
        return True

    def histogram(self):
        return self._sHisto
    
    def buckets(self):
        return [i / self._divids for i in range(-self._half, self._half+1)]


    