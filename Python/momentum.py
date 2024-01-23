import math

class Momentum:
    __slots__ = ('_sFactor', '_maxOrder', '_binding')

    def __init__(self, maxOrder:int=200, binding:int=5) -> None:
        self._maxOrder = maxOrder
        self._binding = binding
        self._sFactor = [0] * maxOrder

        divid = 32
        divid2 = divid * divid
        limit = divid * binding
        for i in range(-limit, limit + 1):
            x2 = i*i / divid2
            pdf = 1.0 / math.sqrt(2*math.pi) * math.exp(- x2 * 0.5) / divid
            sq = 1
            for j in range(maxOrder):
                self._sFactor[j] += pdf * sq
                sq *= x2

    def factor(self, n:int):
        if (n % 2) == 1:
            return 0
        n //= 2
        if n > self._maxOrder:
            return 0
        return self._sFactor[n]



