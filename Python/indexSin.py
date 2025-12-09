import enum
import math
import os

import varDbl


class SinSource (enum.StrEnum):
    Quart = 'Quart',
    Full = 'Full',
    Fixed = 'Fixed',
    Limit = 'Limit',
    Lib = 'Lib',
    Prec = 'Prec'


class IndexSin:
    '''
    A sin() with index frequence as input, with resolution PI /(1 << "order)
    '''
    MIN_ORDER = 1
    MAX_ORDER = 18

    _size = 1 << MAX_ORDER
    _half = 1 << (MAX_ORDER - 1)

    _sSinQuart = None
    _sSinFull = None
    _sSinFixed = None

    __slots__ = ['_order', '_sinSource', '_sSin', '_sCos', '_order']

    @staticmethod
    def validateOrder(order:int):
        if not (IndexSin.MIN_ORDER <= order <= IndexSin.MAX_ORDER):
            raise RuntimeError(f'order={order} is not in the range of [{IndexSin.MIN_ORDER}, {IndexSin.MAX_ORDER}]')

    @staticmethod
    def validateSize(size:int) -> int:
        for order in range(IndexSin.MIN_ORDER, IndexSin.MAX_ORDER + 1):
            if (1 << order) == size:
                return order
        else:
            raise RuntimeError(f'Invalid input array size {size} which is not 2^{order}')


    def __init__(self, sinSource = SinSource.Quart,
                 sCosSin:tuple[varDbl.VarDbl]=None) -> None:                
        self._sinSource = sinSource
        match sinSource:
            case SinSource.Prec:
                IndexSin._sSinPrec = IndexSin.read(SinSource.Prec)
                self._order = IndexSin.MAX_ORDER
                self._sSin = IndexSin._sSinPrec
                self._sCos =  None
            case SinSource.Quart:
                if not IndexSin._sSinQuart:
                    quart = IndexSin._half >> 1
                    sSin = [math.sin(math.pi * i /IndexSin._size) for i in range(quart)] + \
                        [math.cos(math.pi * (quart - i) / IndexSin._size) for i in range(quart + 1)]
                    IndexSin._sSinQuart = tuple([varDbl.VarDbl(v) for v in sSin])
                self._order = IndexSin.MAX_ORDER
                self._sSin = IndexSin._sSinQuart
                self._sCos =  None
            case SinSource.Full:
                if not IndexSin._sSinFull:
                    sSin = [math.sin(math.pi * i /IndexSin._size) for i in range(IndexSin._size)]
                    IndexSin._sSinFull = tuple([varDbl.VarDbl(v) for v in sSin])
                self._order = IndexSin.MAX_ORDER
                self._sSin = IndexSin._sSinFull
                self._sCos =  None
            case SinSource.Fixed:
                if not IndexSin._sSinFixed:
                    sSin = [math.sin(math.pi * i /IndexSin._size) for i in range(IndexSin._size)]
                    IndexSin._sSinFixed = tuple([varDbl.VarDbl(v, math.ulp(1.)) for v in sSin])
                self._order = IndexSin.MAX_ORDER
                self._sSin = IndexSin._sSinFixed
                self._sCos =  None
            case SinSource.Limit:
                if not sCosSin:
                    raise RuntimeError('No sCosSin for SinSource.Limit')
                self._order = IndexSin.validateSize(len(sCosSin) >> 1)
                self._sCos = sCosSin[::2] + [-sCosSin[0]]
                self._sSin = sCosSin[1::2]
            case SinSource.Lib:
                self._order = None
                self._sSin = None
                self._sCos =  None
            case _:
                raise RuntimeError(f'sinSource={sinSource}')

    @staticmethod
    def header():
        return 'Index\tX\tValue\tUncertainty\n'

    def dump(self, order:int=MAX_ORDER) -> None:
        IndexSin.validateOrder(order)
        with open(f'./Python/Output/IndexSin_{self._sinSource}_{order}.txt', 'w') as f:
            f.write(IndexSin.header())
            size = 1 << order
            for i in range(size + 1):
                v = self.sin(i, order)
                f.write(f'{i}\t{i/size:.20e}\t{v.value():.20e}\t{v.uncertainty():.20e}\n')

    @staticmethod
    def read(sinSource:SinSource, order:int=MAX_ORDER, errorFold:float=1e6) -> list[varDbl.VarDbl]:
        '''
        {errorFold} indicates the difference between C++ long double and python double
        '''
        IndexSin.validateOrder(order)
        size = 1 << order
        if os.getcwd().endswith('\\Python'):
            dumpPath = f'../Cpp/Output/IndexSin_{sinSource}_{order}.txt'
        elif os.getcwd().endswith('\\VarianceArithmetic'):
            dumpPath = f'./Cpp/Output/IndexSin_{sinSource}_{order}.txt'
        else:
            raise RuntimeError(f'Invalid current working directory {os.getcwd()}')
        sSin = []
        with open(dumpPath) as f:
            hdr = next(f)
            if hdr != IndexSin.header():
                raise RuntimeError(f'Invalid header in file {dumpPath}: {hdr} vs {IndexSin.header()}')
            for ln, line in enumerate(f):
                sWord = line.split('\t')
                if ln != int(sWord[0]):
                    raise RuntimeError(f'Invalid index {sWord[0]} vs {ln} in {dumpPath}: {line}')
                phase = ln / size
                if float(sWord[1]) != phase:
                    raise RuntimeError(f'Invalid phase {sWord[1]} vs {phase} for index {ln} in {dumpPath}: {line}')
                sin = math.sin(math.pi * phase)
                val, unc = map(float, sWord[2:])
                if unc < 0:
                    raise RuntimeError(f'Invalid uncertainty {sWord[3]} for index {ln} in {dumpPath}: {line}')
                if abs(val - sin) > unc * errorFold:
                    if unc > 0:
                        raise RuntimeError(f'Invalid value {sWord[2]} and uncertainty {sWord[3]} for sin={sin} diffFold={(val - sin)/unc} at index {ln} in {dumpPath}: {line}')
                sSin.append(varDbl.VarDbl(val, unc))
            if len(sSin) < (size >> 1) + 1:
                raise RuntimeError(f'Invalid sin count {len(sSin)} < {(size >> 1) + 1} in {dumpPath}: {line}')
            return sSin

    @property  
    def sinSource(self) -> int:
        return self._sinSource
        
    def get_index(self, freq:int, order:int) ->int:
        '''
        get index into _sSin, with -index means -sin
        '''
        size = 1 << order
        div = freq // size
        rem = freq % size
        if ((self.sinSource == SinSource.Prec) or (self.sinSource == SinSource.Quart)) and (rem > (size >> 1)):
            rem = size - rem
        if div & 1:
            return - rem
        else:       
            return rem


    def sin(self, freq:int, order:int) -> varDbl.VarDbl:
        if self.sinSource == SinSource.Lib:
            v = math.sin(math.pi * freq / (1 << order))
            return varDbl.VarDbl(v)
        elif self.sinSource == SinSource.Limit:
            if self._order < order:
                raise RuntimeError(f'For SinSource.Limit, {self._order} < {order}')
            idx = self.get_index(freq, order)
            if idx >= 0:
                return self._sSin[idx << (self._order - order)]
            else:
                return -self._sSin[(-idx) << (self._order - order)]
        else:
            IndexSin.validateOrder(order)
            idx = self.get_index(freq << (IndexSin.MAX_ORDER - order), IndexSin.MAX_ORDER)
            if idx >= 0:
                return self._sSin[idx]
            else:
                return -self._sSin[-idx]

    def cos(self, freq:int, order:int) -> varDbl.VarDbl:
        if self.sinSource == SinSource.Lib:
            v = math.cos(math.pi * freq / (1 << order))
            return varDbl.VarDbl(v)
        elif self.sinSource == SinSource.Limit:
            if self._order < order:
                raise RuntimeError(f'For SinSource.Limit, {self._order} < {order}')
            '''
            0->7:0->7, 8:8, 9->15:7->1, 16->23:0->7, 24:8    
            -1->-7:1-7, -8:8, -9->-15:7->1, -16->23:0->7, -24:8
            '''
            freq = abs(freq)
            freq <<= self._order - order
            size = 1 << self._order
            div = freq // size
            rem = freq % size
            if div & 1:
                rem = size - rem
            return self._sCos[rem]
        else:
            IndexSin.validateOrder(order)
            return self.sin(freq + (1 << (order - 1)), order)
      
        

