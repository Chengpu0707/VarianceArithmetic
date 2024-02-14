
import abc
import enum
import math

from indexSin import IndexSin



class FFTSinSource (enum.Enum):
    INDEX_SIN = 1,
    LIB_SIN = 2,
    LIMITED_SIN = 3
    UNCERTAIN_SIN = 4


class FFTBase (abc.ABC):
    MAX_ORDER = 18
    _bitReversedIndex = {}

    @abc.abstractmethod
    def source(self) ->FFTSinSource:
        '''
        sin to use for FFT phase
        '''

    @abc.abstractmethod
    def sin(self, freq:int, order:int) -> None:
        '''
        sin(Math.pi * freq / (1 << order))
        '''

    @abc.abstractmethod
    def cos(self, freq:int, order:int) -> None:
        '''
        cos(Math.pi * freq / (1 << order))
        '''

    @staticmethod
    def bitReversedIndices(order:int) -> tuple[int]:
        '''
        Return bit reversed indices
        from NumericalRecipesinC.pdf
        '''
        if sRes := FFTBase._bitReversedIndex.get(order):
            return sRes
        N = 1 << order
        sRes = [0] * N
        M = N >> 1
        j = 0
        for i in range(N):
            sRes[i] = j
            k = M 
            while (k != 0) and (j >= k):
                j -= k
                k >>= 1
            j += k
        FFTBase._bitReversedIndex[order] = sRes
        return sRes
    
    def transform(self, sData:list, forward:bool, useOriginalArray:bool=False):
        '''
        "sData": an array of size (2<<order), with each datum contains (real, image)
        "forward": true for forware transformation, false for backward transformation.
        "useOriginalArray": change sData without allocate a result array
        '''
        for order in range(2, FFTBase.MAX_ORDER):
            if (2 << order) == len(sData):
                break
        if order > FFTBase.MAX_ORDER:
            raise RuntimeError(f'Invalid input array size {len(sData)} which is not 2^{order}')

        sRes = sData if useOriginalArray else [0] * (2 << order)
        
        sIndex = FFTBase.bitReversedIndices(order)
        for i in range(len(sIndex)):
            j = sIndex[i]
            sRes[i << 1], sRes[(i << 1) + 1] = sData[j << 1], sData[(j << 1) + 1]

        for i in range(0, len(sIndex), 2):
            j = i << 1
            sRes[j], sRes[j + 2] = sRes[j] + sRes[j + 2], sRes[j] - sRes[j + 2]
            sRes[j + 1], sRes[j + 3] = sRes[j + 1] + sRes[j + 3], sRes[j + 1] - sRes[j + 3]

        for o in range(2, order + 1):
            k = 1 << o
            for j in range(k >> 1):
                cos = self.cos(j, o)
                sin = self.sin(j, o) if forward else - self.sin(j, o)
                for i in range(0, len(sIndex), k):
                    i0 = (i + j) << 1
                    i1 = i0 + k
                    rd = sRes[i1] * cos - sRes[i1 + 1] * sin
                    id = sRes[i1] * sin + sRes[i1 + 1] * cos
                    sRes[i0], sRes[i1] = sRes[i0] + rd, sRes[i0] - rd
                    sRes[i0 + 1], sRes[i1 + 1] = sRes[i0 + 1] + id, sRes[i0 + 1] - id

        if not forward:
            for i in range(len(sRes)):
                sRes[i] /= len(sIndex)

        return sRes


class FFTIndexSin (FFTBase):
    _indexSin = IndexSin(FFTBase.MAX_ORDER - 1)

    def source(self):
        return FFTSinSource.INDEX_SIN
    
    def sin(self, freq:int, order:int):
        return FFTIndexSin._indexSin.sin(freq *(1 <<(FFTBase.MAX_ORDER - order)))
    
    def cos(self, freq:int, order:int):
        return FFTIndexSin._indexSin.cos(freq *(1 <<(FFTBase.MAX_ORDER - order)))


class FFTLibSin (FFTBase):
    def source(self):
        return FFTSinSource.INDEX_SIN
    
    def sin(self, freq:int, order:int):
        return math.sin(math.pi *freq /(1 << (order - 1)))
    
    def cos(self, freq:int, order:int):
        return math.cos(math.pi *freq /(1 << (order - 1)))
    

class FFTLimitedSin (FFTBase):
    _indexSin = IndexSin(FFTBase.MAX_ORDER - 1)

    def source(self):
        return FFTSinSource.LIMITED_SIN
    
    def sin(self, freq:int, order:int):
        idx = FFTLimitedSin._indexSin.get_index(freq)
        return math.sin(math.pi *idx /(1 << (order - 1))) if idx >= 0 else\
              -math.sin(math.pi *-idx /(1 << (order - 1)))
    
    def cos(self, freq:int, order:int):
        idx = FFTLimitedSin._indexSin.get_index(freq)
        return math.cos(math.pi *idx /(1 << (order - 1))) if idx >= 0 else\
              -math.sin(math.pi *idx /(1 << (order - 1)))


class FFTUncertainSin (FFTBase):
    _indexSin = None

    def __init__(self) -> None:
        super().__init__()
        if not FFTUncertainSin._indexSin:
            FFTUncertainSin._indexSin = IndexSin(FFTBase.MAX_ORDER - 1)
            err = FFTUncertainSin._indexSin.withUncertainty()
            if err:
                raise ValueError(err)

    def source(self):
        return FFTSinSource.UNCERTAIN_SIN
    
    def sin(self, freq:int, order:int):
        return FFTUncertainSin._indexSin.sin(freq *(1 <<(FFTBase.MAX_ORDER - order)))
    
    def cos(self, freq:int, order:int):
        return FFTUncertainSin._indexSin.cos(freq *(1 <<(FFTBase.MAX_ORDER - order)))
