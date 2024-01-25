import math

from varDbl import VarDbl

class IndexSin:
    '''
    A sin() with index frequence as input, with resolution \pi/size()
    size() = 1 << "order".
    When the index "freq" == size, the result is sin(\pi).
    This function has minimal float error when "withUncertainty"==False.

    When "withUncertainty"==True, use regression to calculate the sin
    '''
    
    __slots__ = ['_size', '_half', '_sSin']

    def __init__(self, order=18, withUncertainty:bool=False) -> None:
        if order < 2:
            raise ValueError(f'order {order} is less than 4 for IndexSin')
        self._size = 1 << order
        self._half = self._size >> 1
        if withUncertainty:
            raise NotImplemented() 
        else:
            self._sSin = [0] + [math.sin(i/self._size *math.pi) for i in range(1, self._half >> 1)] + \
                    [math.cos((1/4 - i/self._size)*math.pi) for i in range(self._half >> 1)] + [1] 
        
    def size(self):
        return self._size
        
    def _get_sin_index(self, freq:int) ->int:
        '''
        get index into _sSin, with -index means -sin
        '''
        div = freq // self._half
        rem = freq % self._half
        if div & 1:
            div += 1
            rem -= self._half
        if div & 2:
            return -rem
        else:       
            return rem

    def _get_cos_index(self, freq:int) ->int:
        return self._get_sin_index(freq + self._half)
    
    def sin(self, freq:int) -> float:
        idx = self._get_sin_index(freq)
        return self._sSin[idx] if idx >= 0 else -self._sSin[-idx]

    def cos(self, freq:int) -> float:
        idx = self._get_cos_index(freq)
        return self._sSin[idx] if idx >= 0 else -self._sSin[-idx]

    def tan(self, freq:int) -> float:
        return self.sin(freq) / self.cos(freq)


    


