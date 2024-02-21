
import math

class Stat:
    __slots__ = ['_count', '_min', '_max', '_sum', '_sum2']

    def __init__(self) -> None:
        self._count = 0
        self._min = 0.0
        self._max = 0.0
        self._sum = 0.0
        self._sum2 = 0.0

    def accum(self, value:float) ->bool:
        value = float(value)
        if not math.isfinite(value):
            return False
        if (not self._count) or (self._min > value):
            self._min = value
        if (not self._count) or (self._max < value):
            self._max = value
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

    def range(self):
        return self._devs
    
    def less(self):
        return self._less

    def more(self):
        return self._more

    def stat(self):
        return self._stat
    
    def accum(self, value:float) ->bool:
        if not self._stat.accum(value):
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


    