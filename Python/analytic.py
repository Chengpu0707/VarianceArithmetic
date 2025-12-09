import functools
import math
import operator
import sympy
import typing

import momentum

class momentum (sympy.Function):
    '''
    Variance momentum for momentum.Normal distribution
    '''
    is_negative = False
    mmt = momentum.Normal()

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
        #return sympy.Integer(math.prod(range(n - 1, 0, -2)))
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
            maxOrder:int=448) \
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
        sX: tuple[sympy.Symbol], 
        maxOrder:int=448) \
    -> tuple[dict[typing.Union[tuple[int],int], sympy.Function]]:
    '''
    Calculate Taylor expansion of {func} to a set of variables {sX}, 
        up to {maxOrder} of Taylor expansion.
    return the Taylor series, the bias series, and the variance series,
        within each indexed as the power of the variable set
    '''
    if not sX:
        raise ValueError(f'No variable {sX} for {func}')
    if not _is_iterable(sX):
        raise ValueError(f'Variable {sX} is not iterable for {func}')
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
            if sum(key) >= maxOrder:
                continue
            if key not in sVar:
                sVar[key] = 0
            sVar[key] += sDiff[key1] * sDiff[key2] * _key_momentum(key)
            if [k for k in key1 if (k % 2)] and [k for k in key2 if (k % 2)]:
                continue
            sVar[key] -= sDiff[key1] * _key_momentum(key1) * sDiff[key2] * _key_momentum(key2)
    sDiff[tuple([0]*len(sX))] = func

    return sDiff, sBias, sVar


def taylor(func: sympy.Function, 
        sXnVar: tuple[tuple[sympy.Symbol]], 
        maxOrder:int=448,
        minTermimationOrder:int=20) \
    -> tuple[sympy.Function]:
    '''
    Calculate Taylor expansion of {func} to a set of variables {sXnVar}, 
        up to {maxOrder} of Taylor expansion
    The input is a tupel for mutiple dimension of input.
    Each input contains:
        symbolic x, 
        variance(x), which could either be a symbol or a value
        optional value(x).
    The result is four part: sDiff, bias, variance, convergency, 
        with convergency be either true or false or None
    '''
    if not sXnVar:
        raise ValueError(f'No variable {sXnVar} for {func}')
    if _is_iterable(sXnVar):
        sX = []
        sVarX = []
        sValX = []
        for i, XnVar in enumerate(sXnVar):
            if not _is_iterable(sXnVar):
                raise ValueError(f'Invalid #{i} variable {XnVar} in {sXnVar} for {func}')
            match len(XnVar):
                case 0:
                    sX.append(XnVar)
                    sVarX.append(0)
                    sValX.append(None)
                case 1:
                    sX.append(XnVar[0])
                    sVarX.append(0)
                    sValX.append(None)
                case 2:
                    sX.append(XnVar[0])
                    sVarX.append(XnVar[1])
                    sValX.append(None)
                case 3:
                    sX.append(XnVar[0])
                    sVarX.append(XnVar[1])
                    sValX.append(XnVar[2])
                case _:
                    raise ValueError(f'Invalid # {i} input variable {sXnVar} for {func}')
    else:
        sX = (sXnVar,)
        sVarX = (0,)
        sValX = (None,)
    assert len(sX) == len(sVarX) == len(sValX)
    sDiff, sBias, sVar = taylor_series(func, sX, maxOrder=maxOrder)
    sSub = {x: val for x, val in zip(sX, sValX) if val is not None}
    if sSub:
        sDiff = {k: v.subs(sSub) for k,v in sDiff.items()}
        sBias = {k: v.subs(sSub) for k,v in sBias.items()}
        sVar =  {k: v.subs(sSub) for k,v in sVar.items()}
    bias = 0
    for k,v in sBias.items():
        bias += v.evalf() * functools.reduce(operator.mul, [x**pw for pw,x in zip(k, sVarX)])
    sExpand = [0]
    for k,v in sVar.items():
        idx = sum(k) //2
        if len(sExpand) <= idx:
            sExpand.extend([0] * (idx - len(sExpand) + 1))
        sExpand[idx] += v.evalf() * functools.reduce(operator.mul, [x**(pw//2) for pw,x in zip(k, sVarX)])
    convergency = True
    if (minTermimationOrder > 0) and (len(sDiff) >= momentum.Normal.MAX_ORDER):
        for i in range(1, minTermimationOrder //2):
            try:
                if bool((sExpand[-i] / sExpand[-i - 1]) > 1):
                    convergency = False
                    break
            except TypeError:
                convergency = None
                break
    return sDiff, bias, sum(sExpand), convergency






