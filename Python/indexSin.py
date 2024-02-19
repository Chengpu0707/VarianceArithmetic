import bisect
import math
import os
import typing

from varDbl import VarDbl

class IndexSin:
    '''
    A sin() with index frequence as input, with resolution \pi/size()
    size() = 1 << "order".
    When the index "freq" == size, the result is sin(\pi).
    This function has minimal float error when "withUncertainty"==False.

    When "withUncertainty"==True, use regression to calculate the sin
    '''
    zero = VarDbl(0,0)
    one = VarDbl(1,0)
    two = VarDbl(2,0)
    half = VarDbl(1/2,0)

    staticmethod
    def path(order:int) ->str:
        return f'./Python/Output/UncertainSin_{order}.txt'
    
    header = "Order\tIndex"\
             "\tSin Value\tSin Uncertainty\tSin Normalized Error"\
             "\tCos Value\tCos Uncertainty\tCos Normalized Error"\
             "\tError Value\tError Uncertainty\tNormalized Error\n"
    
    __slots__ = ['_order', '_size', '_half', '_sSin']

    def __init__(self, order) -> None:
        if order < 3:
            raise ValueError(f'order {order} is less than 4 for IndexSin')
        self._order = order
        self._size = 1 << order
        self._half = self._size >> 1
        self._sSin = [0] + [math.sin(i/self._size *math.pi) for i in range(1, self._half >> 1)] + \
                [math.cos((1/4 - i/self._size)*math.pi) for i in range(self._half >> 1)] + [1] 
        
    def withUncertainty(self, incomplete:bool=False) -> typing.Optional[str]:
        '''
        Return None for successfully reading the file which contain sin with uncertainty.
        Otherwise return error string and restore the original float/int values

        "binding" is for comparing the z-difference between lib sin/cos and the calculated sin/cos
        '''
        sSin = self._sSin[:]
        with open (IndexSin.path(self._order)) as f:
            title = next(f)
            if title != IndexSin.header:
                return f'Wrong header: {title}'
            for ln, line in enumerate(f):
                sWord = line.strip().split('\t')
                if len(sWord) != 11:
                    self._sSin = sSin
                    return f'line #{ln} is invalid: {line}'
                ord, idx = map(int, sWord[:2])
                if not (0 <= ord <= self._order):
                    self._sSin = sSin
                    return f'line #{ln} order={ord} outside the range of [0, {self._order}]'
                if not (0 <= idx <= self._half):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} outside the range of [0, {self._half}]'
                arc = math.pi * idx/self._size
                sValue = list(map(float, sWord[2:])) 
                if sValue[2] if (sValue[1] <= 0) else \
                        math.ulp(sValue[2]) < abs((sValue[0] - math.sin(arc)) / sValue[1] - sValue[2]):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} sin values disagree: {math.sin(arc)} vs {sValue[0]}~{sValue[1]} = {sValue[2]}'
                if sValue[4] if (sValue[4] <= 0) else \
                        math.ulp(sValue[3]) < abs((sValue[3] - math.cos(arc)) / sValue[4] - sValue[5]):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} cos values disagree: {math.cos(arc)} vs {sValue[3]}~{sValue[4]} = {sValue[5]}'
                if sValue[7] if (sValue[7] <= 0) else \
                        math.ulp(sValue[8]) < abs(sValue[6] / sValue[7] - sValue[8]):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} error values disagree: {sValue[8]} vs {sValue[6]}~{sValue[7]}'
                self._sSin[idx] = VarDbl(sValue[0], sValue[1])
                self._sSin[self._half - idx] = VarDbl(sValue[3], sValue[4])
                error = VarDbl(sValue[6], sValue[7])
                err = self._sSin[idx] * self._sSin[idx] + \
                        self._sSin[self._half - idx] * self._sSin[self._half - idx] - 1
                if error.uncertainty() < math.ulp(sValue[6]) < abs(error.value() - err.value()):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} error values disagree: {error} vs {err}'
                '''
                if sValue[6] if (sValue[6] <= 0) else \
                        math.ulp(sValue[6]) < abs(error.value() - err.value()):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} error values disagree: {error} vs {err}'
                if sValue[7] if (sValue[7] <= 0) else \
                        math.ulp(sValue[7]) < abs(error.uncertainty() - err.uncertainty()):
                    self._sSin = sSin
                    return f'line #{ln} index={idx} error uncertainties disagree: {error} vs {err}'
                '''
        sMissing = [i for i, e in enumerate(self._sSin) if not isinstance(e, VarDbl)]
        if sMissing and (not incomplete):
            self._sSin = sSin
            return f'order {self._order} missing {len(sMissing)} indices: {sMissing}'

        return None

    def calc(self):
        '''
        Return None for successfully calculating the file which contain sin with uncertainty.
        Otherwise return error string

        "binding" is for comparing the z-difference between lib sin/cos and the calculated sin/cos
        '''
        try:
            self.withUncertainty(incomplete=True)
        except BaseException as ex:
            print(f'Fail to read {IndexSin.path(self._order)}: {ex}')
        self._sSin[0] = IndexSin.zero
        self._sSin[self._half] = IndexSin.one
        exist = os.path.isfile(IndexSin.path(self._order))
        with open (IndexSin.path(self._order), 'a' if exist else 'w') as f:
            if not exist:
                f.write(IndexSin.header)
                f.write('0\t0\t0\t0\t0\t1\t0\t0\t0\t0\t0\n')
            self._calc(0, self._half, 1, f)

    def _calc(self, begin, end, order, file):
        if order >= self._order:
            return None
        smid = (begin + end) >> 1
        cmid = self._half - smid
        if (type(self._sSin[smid]) != VarDbl) or (type(self._sSin[cmid]) != VarDbl):
            x = self._sSin[self._half - begin] * self._sSin[self._half - end] \
                - self._sSin[begin] * self._sSin[end]
            self._sSin[smid] = ((IndexSin.one - x) *IndexSin.half) ** 0.5
            self._sSin[cmid] = ((IndexSin.one + x) *IndexSin.half) ** 0.5
            arc = math.pi * smid/self._size
            err = self._sSin[smid] **2 + self._sSin[cmid] **2 - 1   
            file.write(f'{order}\t{smid}')
            file.write(f'\t{self._sSin[smid].value():.16e}\t{self._sSin[smid].uncertainty():.16e}\t{(self._sSin[smid].value() - math.sin(arc))/self._sSin[smid].uncertainty():.16e}')
            file.write(f'\t{self._sSin[cmid].value():.16e}\t{self._sSin[cmid].uncertainty():.16e}\t{(self._sSin[cmid].value() - math.cos(arc))/self._sSin[cmid].uncertainty():.16e}')
            file.write(f'\t{err.value():.16e}\t{err.uncertainty():.16e}\t{err.value()/err.uncertainty():.16e}\n')
            file.flush()
        self._calc(begin, smid, order + 1, file)
        self._calc(smid, end, order + 1, file)
        
    def size(self):
        return self._size
        
    def get_index(self, freq:int) ->int:
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

    def sin(self, freq:int) -> float:
        idx = self.get_index(freq)
        return self._sSin[idx] if idx >= 0 else -self._sSin[-idx]

    def cos(self, freq:int) -> float:
        return self.sin(freq + self._half)

    def tan(self, freq:int) -> float:
        return self.sin(freq) / self.cos(freq)
      
    def arc_sin(self, value:float) -> float:
        if not -1 <= value <= 1:
            raise ValueError(f'Invalid sine value {value}')
        sign = 1 if value >=0 else -1
        value = abs(value)
        idx = bisect.bisect_left(self._sSin, value)
        if idx == len(self._sSin) - 1:
            return sign * idx 
        return sign * (idx + (value - self._sSin[idx]) / (self._sSin[idx + 1] - self._sSin[idx]))
        



    


