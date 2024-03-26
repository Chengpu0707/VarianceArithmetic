import fractions
import math
import itertools
import random
import typing

import varDbl

ElementTypes = (int, float, fractions.Fraction, varDbl.VarDbl)
#TODO: ElementType = typing.Union[*ElementTypes] # Unpack is not allowed in this context
ElementType = typing.Union[int, float, fractions.Fraction, varDbl.VarDbl]




def permutSign(sPermut:tuple[int]) -> int:
    '''
    According to https://en.wikipedia.org/wiki/Parity_of_a_permutation,
    return +1 or -1.

    Other more efficient way: https://stackoverflow.com/questions/20702782/efficiently-determine-the-parity-of-a-permutation
    '''
    cnt = 0
    for i in range(len(sPermut) - 1):
        for j in range(i + 1, len(sPermut)):
            if sPermut[i] > sPermut[j]:
                cnt += 1
    return -1 if (cnt % 2) else 1


def isSquareMatrix(ssMatrix:tuple[tuple[ElementType]], sType=ElementTypes) ->bool:
    '''
    Not using numpy for a matrix to avoid coercing int to int32.
    Use tuple to avoid changing of the matrix size.
    "sType" holds the allowed element type
    '''
    if type(ssMatrix) != tuple:
        return False
    size = len(ssMatrix)
    if [1 for sMatrix in ssMatrix if (len(sMatrix) != size) or type(sMatrix) != tuple]:
        return False
    if [1 for sMatrix in ssMatrix for val in sMatrix if type(val) not in sType]:
        return False
    return True


def createIntMatrix(size:int, randRange) ->tuple[tuple[int]]:
    '''
    Create an int matrix of size "size" with each element uniformly distributed between [-"randRange", +"randRange"]
    '''
    size = abs(size)
    return tuple([tuple([random.randint(-randRange, +randRange) for col in range(size)]) 
                  for row in range(size)])

def createHilbertMatrix(size:int) ->tuple[tuple[fractions.Fraction]]:
    '''
    https://en.wikipedia.org/wiki/Hilbert_matrix
    '''
    size = abs(size)
    return tuple([tuple([fractions.Fraction(1, row + col + 1) for col in range(size)]) 
                  for row in range(size)])

def addNoise(ssMatrix:tuple[tuple[typing.Union[int]]], noise:float) -> tuple[tuple[varDbl.VarDbl]]:
    '''
    Add Gaussian "noise" to "ssMatrix", with the actual noise is adjusted by ELEMENT_RANGE.
    '''
    if not isSquareMatrix(ssMatrix, sType=(int,fractions.Fraction)):
        raise ValueError(f'Invalid int or float')
    size = len(ssMatrix)
    return tuple([tuple([varDbl.VarDbl(ssMatrix[row][col] + noise * random.normalvariate(), 
                                       noise**2 + math.ulp(ssMatrix[row][col])**2, True) 
                         for col in range(size)]) 
                  for row in range(size)])


def linear(ssMatrix:tuple[tuple[ElementType]], scale:ElementType=1, offset:ElementType=0) -> tuple[tuple[ElementType]]:
    if not isSquareMatrix(ssMatrix):
        raise ValueError(f'The input square matrix is illegal for linear(): {ssMatrix}')
    if type(scale) not in ElementTypes:
        raise ValueError(f'The input scale {scale} is illegal type {type(scale)} for linear()')
    size = len(ssMatrix)
    return tuple([tuple([ssMatrix[i][j] * scale + offset for j in range(size)]) 
                  for i in range(size)])
            

def multiply(ssMatrix1:tuple[tuple[ElementType]], ssMatrix2:tuple[tuple[ElementType]]) -> tuple[tuple[ElementType]]:
    if not isSquareMatrix(ssMatrix1):
        raise ValueError(f'The input square matrix 1 is illegal for multiply(): {ssMatrix1}')
    if not isSquareMatrix(ssMatrix2):
        raise ValueError(f'The input square matrix 1 is illegal for multiply(): {ssMatrix2}')
    size = len(ssMatrix1)
    if size != len(ssMatrix2):
        raise ValueError(f'The input square matrix 1 and 2 has different size {size} vs {len(ssMatrix2)}')
    return tuple([tuple([sum([ssMatrix1[i][k] * ssMatrix2[k][j] for k in range(size)]) for j in range(size)]) 
                  for i in range(size)])
            

def adjugate(ssMatrix:tuple[tuple[ElementType]]) -> tuple[ElementType, tuple[tuple[ElementType]]]:
    '''
    Calculate determinant and the adjugate matrix for "ssMatrix".
    If the ElementType contains mixed types, the result promotion is int -> Fraction -> float.
    '''
    if not isSquareMatrix(ssMatrix):
        raise ValueError(f'The input square matrix is illegal for determinant(): {ssMatrix}')
    size = len(ssMatrix)
    if size == 1:
        return ssMatrix[0][0], ssMatrix
    elif size == 2:
        return ssMatrix[0][0]*ssMatrix[1][1] - ssMatrix[0][1]*ssMatrix[1][0], \
               tuple([tuple([ssMatrix[1][1], -ssMatrix[0][1]]), tuple([-ssMatrix[1][0], ssMatrix[0][0]])])

    sPermut = {permut: permutSign(permut) for permut in itertools.permutations(range(size), size)}

    value = 0
    variance = 0
    sCofVar = {(i,j):0 for i in range(size) for j in range(size)}
    for permut, sign in sPermut.items():
        val = sign
        var = 1
        for x, y in enumerate(permut):
            if type(ssMatrix[x][y]) == varDbl.VarDbl:
                val *= ssMatrix[x][y].value()
                var *= ssMatrix[x][y].variance()
            else:
                val *= ssMatrix[x][y]
                var = 0
        value += val
        variance += var
        for x, y in sCofVar:
            if (permut[x] == y) and (var > 0):
                sCofVar[(x,y)] += var / ssMatrix[x][y].variance()
    
    sCofVal = {}
    for m in range(1, size):
        sVal = {}
        sVar = {}
        for sX in itertools.combinations(range(size), m):
            for permut, sign in sPermut.items():
                val = sign
                var = 1
                sY = tuple([(x,permut[x]) for x in sX])
                for x, y in enumerate(permut):
                    if x in sX:
                        var *= ssMatrix[x][y].variance() if type(ssMatrix[x][y]) == varDbl.VarDbl else 0
                    else:
                        val *= ssMatrix[x][y].value() if type(ssMatrix[x][y]) == varDbl.VarDbl else ssMatrix[x][y]
                if sY in sVal:
                    sVal[sY] += val
                    assert sVar[sY] == var
                else:
                    sVal[sY] = val
                    sVar[sY] = var
        for sY, var in sVar.items():
            try:
                var *= sVal[sY] ** 2
                variance += var
            except OverflowError as ex:
                print(f'adjugate(): size={size}, m={m}, max_element={max([max(abs(val.value()) for val in row) for row in ssMatrix])}')
                raise ex
            if var > 0:
                for x, y in sY:
                    sCofVar[(x, y)] += var / ssMatrix[x][y].variance()
        if m == 1:
            for k, v in sVal.items():
                sCofVal[k[0]] = v
                sCofVar[k[0]] -= v**2
 
    return varDbl.VarDbl(value, variance, True) if variance > 0 else value, \
           tuple([tuple([varDbl.VarDbl(sCofVal[(i,j)], sCofVar[(i,j)], True) if 0 < sCofVar[(i,j)] else sCofVal[(i,j)] 
                         for i in range(size)]) 
                  for j in range(size)])



