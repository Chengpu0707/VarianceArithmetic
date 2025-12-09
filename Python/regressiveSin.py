'''
Regressive calculation of sine and cosine from sin(0) = 0 and sin(PI/2) = 1
'''
import os

from indexSin import IndexSin, SinSource
import varDbl


class RegressiveSin:
    '''
    Use regression to generate sin(j /(1<<order) *math.pi)
    '''
    ZERO = varDbl.VarDbl(0,0)
    ONE = varDbl.VarDbl(1,0)
    HALF = varDbl.VarDbl(1/2,0)

   
    def __init__(self, order) -> None:
        IndexSin.validateOrder(order)
        self._order = order
        self._size = 1 << self._order
        self._half = self._size >> 1
        self._sSin = [None] * (self._half + 1)
        self._indexSin = IndexSin(SinSource.Quart)

    def sin(self, freq:int) -> varDbl.VarDbl:
        return self._indexSin.sin(freq, self._order)

    def calc(self):
        '''
        Return None for successfully calculating the fsin which contain sin with uncertainty.
        Otherwise return error string
        '''
        self._sSin[0] = RegressiveSin.ZERO
        self._sSin[self._half] = RegressiveSin.ONE
        if os.getcwd().endswith('/VarianceArithmetic'):
            dirPath = './Python/Output'
        elif os.getcwd().endswith('/VarianceArithmetic/Python'):
            dirPath = './Output'
        elif os.getcwd().endswith('\\VarianceArithmetic'):
            dirPath = './Python/Output'
        elif os.getcwd().endswith('\\VarianceArithmetic\\Python'):
            dirPath = './Output'
        else:
            raise RuntimeError(f'Unexpected cwd {os.getcwd()}')
        sinPath = f'{dirPath}/RegrSin_{self._order}.txt'
        errPath = f'{dirPath}/RegrSin_{self._order}.err.txt'
        with open (sinPath, 'w') as fsin, open (errPath, 'w') as ferr:
            fsin.write("Order\tIndex\tQuart Value\tQuart Uncertainty\tRegr Value\tRegr Uncertainty\tRegr Error\tRegr Normalized\n")
            fsin.write('0\t0\t0\t0\t0\t0\t0\t0\n'
                        f'0\t{self._half}\t1\t0\t1\t0\t0\t0\n')
            ferr.write("Order\tSin Index\tCos Index\tRegr Error\tRegr Uncertainty\tRegr Normalized\tQuart Error\tQuart Uncertainty\tQuart Normalized\n")
            ferr.write(f'0\t0\t{self._half}\t0\t0\t\t0\t0\t\n')
            self._calc(0, self._half, 1, fsin, ferr)
        sMissing = [i for i,s in enumerate(self._sSin) if s is None]
        if sMissing:
            raise RuntimeError(f'order {self._order} missing {len(sMissing)} indices: {sMissing}')


    def _calc(self, begin, end, order, fsin, ferr):
        if order > self._order:
            return None
        smid = (begin + end) >> 1
        cmid = self._half - smid
        if (self._sSin[smid] is None) or (self._sSin[cmid] is None):
            x = self._sSin[self._half - begin] * self._sSin[self._half - end] - self._sSin[begin] * self._sSin[end]
            if self._sSin[smid] is None:
                self._sSin[smid] = ((RegressiveSin.ONE - x) *RegressiveSin.HALF) ** 0.5
            if self._sSin[cmid] is None:
                self._sSin[cmid] = ((RegressiveSin.ONE + x) *RegressiveSin.HALF) ** 0.5
            err = self._sSin[smid].value() - self.sin(smid).value()
            fsin.write(f'{order}\t{smid}\t{self.sin(smid).value()}\t{self.sin(smid).uncertainty()}\t{self._sSin[smid].value()}\t{self._sSin[smid].uncertainty()}\t{err}')
            if self._sSin[smid].uncertainty():
                fsin.write(f'\t{err/self._sSin[smid].uncertainty()}\n')
            else:
                fsin.write(f'\t\n')
            if smid != cmid:
                err = self._sSin[cmid].value() - self.sin(cmid).value()
                fsin.write(f'{order}\t{cmid}\t{self.sin(cmid).value()}\t{self.sin(cmid).uncertainty()}\t{self._sSin[cmid].value()}\t{self._sSin[cmid].uncertainty()}\t{err}')
                if self._sSin[cmid].uncertainty():
                    fsin.write(f'\t{err/self._sSin[cmid].uncertainty()}\n')
                else:
                    fsin.write(f'\t\n')
            fsin.flush()
            ferr.write(f'{order}\t{smid}\t{cmid}')
            err = self._sSin[smid] **2 + self._sSin[cmid] **2 - 1   
            ferr.write(f'\t{err.value()}\t{err.uncertainty()}')
            if err.uncertainty():
                ferr.write(f'\t{err.value()/err.uncertainty()}')
            else:
                ferr.write('\t')
            err = self.sin(smid) **2 + self.sin(cmid) **2 - 1   
            ferr.write(f'\t{err.value()}\t{err.uncertainty()}')
            if err.uncertainty():
                ferr.write(f'\t{err.value()/err.uncertainty()}\n')
            else:
                ferr.write('\t\n')
        self._calc(begin, smid, order + 1, fsin, ferr)
        self._calc(smid, end, order + 1, fsin, ferr)
        

    


