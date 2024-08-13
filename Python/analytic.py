import functools
import operator
import sympy
import typing

from momentum import Momentum

class momentum (sympy.Function):
    '''
    Variance momentum for Normal distribution
    '''
    is_negative = False
    mmt = Momentum()

    @classmethod
    def eval(cls, n):
        if not isinstance(n, sympy.Integer):
            raise TypeError(f'Invalid type {type(n)} for input n={n}')
        if n.is_negative:
            raise TypeError(f'Invalid int value input n={n}')
        if n.is_odd:
            return 0
        
    def _eval_evalf(self, prec):
        '''
        TODO: not find dps argment for evalf() to convert prec
        '''
        n = self.args[0]
        return sympy.Float(momentum.mmt[int(n)])



def _key_momentum(key:tuple[int]) -> sympy.Function:
    return functools.reduce(operator.mul, [momentum(k) for k in key if k > 0])

def _is_iterable(obj) -> bool:
    try:
        iter(obj)
        return True
    except TypeError:
        return False
    

def _taylor_matrix_series(matx: sympy.Matrix, 
            sX: tuple[sympy.Symbol], 
            maxOrder:int=Momentum.MAX_ORDER) \
    -> sympy.Matrix:
    '''
    Calculate Taylor expansion of {func} to a set of variables in the matrix {sX}, 
        up tp {maxOrder} of Taylor expansion
    return the matrix Taylor series, the bias series, and the variance series,
        with each indexed as the power of the variable set
    '''
    f = functools.partial(taylor_series, sX=sX, maxOrder=maxOrder) 
    sMat = matx.applyfunc(f)
    return sMat.applyfunc(lambda e: e[0]), sMat.applyfunc(lambda e: e[1]), sMat.applyfunc(lambda e: e[2])


def taylor_series(func: sympy.Function, 
        sX: typing.Union[tuple[sympy.Symbol], sympy.Symbol], 
        maxOrder:int=Momentum.MAX_ORDER,
        indexInsteadOfTupleFor1d=True) \
    -> tuple[dict[typing.Union[tuple[int],int], sympy.Function]]:
    '''
    Calculate Taylor expansion of {func} to a set of variables {sX}, 
        up to {maxOrder} of Taylor expansion
    return the Taylor series, the bias series, and the variance series,
        withn each indexed as the power of the variable set
    For 1d Taylor expansion, if {indexInsteadOfTupleFor1d} is True:
        *) {sX} can be a symbol instead of a tuple
        *) The return keys are int rather than (int,)
    '''
    if not sX:
        raise ValueError(f'No variable {sX} for {func}')
    if not _is_iterable(sX):
        sX = (sX,)
    if isinstance(func, sympy.Matrix):
        return _taylor_matrix_series(func, sX, maxOrder=maxOrder) 
    sDiff = {tuple([0]*len(sX)): func}
    sBias = {}
    sVar = {}
    for n in range(1, maxOrder):
        sKey = [k for k in sDiff if sum(k) == (n - 1)]
        if not sKey:
            break
        for k in sKey:
            for i, x in enumerate(sX):
                diff = sympy.diff(sDiff[k], x)
                if not diff:
                    continue
                key = list(k)
                key[i] += 1
                key = tuple(key)
                sDiff[key] = diff / key[i]
    
    del sDiff[tuple([0]*len(sX))]
    for key1 in sDiff:
        if not [k for k in key1 if (k % 2)]:
            sBias[key1] = sDiff[key1] * momentum(sum(key1))
        for key2 in sDiff:
            key = tuple(map(operator.add, key1, key2))
            if [k for k in key if (k % 2)]:
                continue
            if key not in sVar:
                sVar[key] = sympy.Integer(0)
            sVar[key] += sDiff[key1] * sDiff[key2] * _key_momentum(key)
            if [k for k in key1 if (k % 2)] and [k for k in key2 if (k % 2)]:
                continue
            sVar[key] -= sDiff[key1] * _key_momentum(key1) * sDiff[key2] * _key_momentum(key2)
    sDiff[tuple([0]*len(sX))] = func

    if (len(sX) == 1) and indexInsteadOfTupleFor1d:
        sDiff = {k[0]: v for k,v in sDiff.items()}
        sBias = {k[0]: v for k,v in sBias.items()}
        sVar = {k[0]: v for k,v in sVar.items()}
    return sDiff, sBias, sVar


def taylor(func: sympy.Function, 
        sXnVar: typing.Union[tuple[tuple[sympy.Symbol]], tuple[sympy.Symbol]], 
        maxOrder:int=126,
        checkConvergence=False,
        checkRelaibility=False) \
    -> tuple[sympy.Function]:
    '''
    Calculate Taylor expansion of {func} to a set of variables {sXnVar}, 
        up to {maxOrder} of Taylor expansion
    Each input is a tuple containing value(x) +/- variance(x).
    The input is a tupel for mutiple dimension of input.
    The result is two value: bias and variance
    If {checkConvergence} is true, make sure that the vairance is monotonic with order
    If {checkRelaibility} is true, make sure that the vairance has precision less than 1/sigma
    '''
    if not sX:
        raise ValueError(f'No variable {sX} for {func}')
    if not _is_iterable(sX):
        sX = (sX,)
        sVarX = (0,)
    if len(sXnVar) == 2:
        sX = (sXnVar[0],)
        sVarX = (sXnVar[1],)
    else:
        sX = []
        sVarX = []
        for i, XnVar in enumerate(sXnVar):
            if not _is_iterable(XnVar):
                sX.append(XnVar)
                sVarX.append(0)
                continue
            match len(XnVar):
                case 1:
                    sX.append(XnVar[0])
                    sVarX.append(0)
                case 2:
                    sX.append(XnVar[0])
                    sVarX.append(XnVar[1])
                case _:
                    raise ValueError(f'Invalid # {i} input variable {sXnVar} for {func}')
    assert len(sX) == len(sVarX)
    sDiff, sBias, sVar = taylor_series(func, sX, maxOrder=maxOrder, indexInsteadOfTupleFor1d=False)
    bias = 0
    for k,v in sBias.items():
        bias += v * functools.reduce(operator.mul, [x**pw for pw,x in zip(k, sVarX)])
    var = 0
    if checkConvergence or checkRelaibility:
        raise NotImplementedError(f'checkConvergence or checkRelaibility are not implemented for Taylor expansion in general')
    for k,v in sVar.items():
        var += v * functools.reduce(operator.mul, [x**pw for pw,x in zip(k, sVarX)])
    return sDiff, bias, var






