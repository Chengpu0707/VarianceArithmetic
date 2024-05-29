import math

import varDbl

class Momentum:
    '''
    Calculate variance momentum.

    Value becomes inf when MAX_ORDER=127
    Variance becomes inf when MAX_ORDER=123
    '''
    BINDING_FACTOR:float = 5.0
    MAX_ORDER:int = 122

    __slots__ = ('_sFactor', '_maxOrder', '_binding')

    def __init__(self, binding:float=BINDING_FACTOR, maxOrder:int=MAX_ORDER) -> None:
        '''
        "binding" is default to 5 for leakage smaller than 5e-7
        "maxOrder" is default to the max 126 before variance overflows

        This call took about 1.3 second
        '''
        self._binding = binding
        self._maxOrder = maxOrder
        self._sFactor = []
        for j in range(maxOrder):
            self._sFactor.append(varDbl.VarDbl())

        divid = 32
        divid2 = divid * divid
        limit = int(math.ceil(divid * binding))
        for i in range(-limit, limit + 1):
            x2 = varDbl.VarDbl(i*i / divid2, 0)
            pdf = varDbl.VarDbl( 1.0 / math.sqrt(2*math.pi)) * varDbl.VarDbl( math.exp(- x2.value() * 0.5) / divid )
            sq = varDbl.VarDbl( 1, 0 )
            for j in range(maxOrder):
                self._sFactor[j] += pdf * sq
                sq *= x2
            # print(i, pdf, self._sFactor[0])

    def factor(self, n:int) -> varDbl.VarDbl:
        '''
        Fetch the variance momentum at "n".

        To be tested by IPyNb/Momentum.ipynb.
        '''
        if (n % 2) == 1:
            return varDbl.VarDbl( 0, 0 )
        n //= 2
        if n > self._maxOrder:
            return IndexError()
        return self._sFactor[n]



