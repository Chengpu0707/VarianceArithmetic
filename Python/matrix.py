import fractions
import itertools
import math
import random
import typing

import histo
import varDbl

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


def isSquareMatrix(ssMatrix:tuple[tuple[ElementType]], 
            sType=(varDbl.VarDbl, fractions.Fraction, int, float)) ->bool:
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


ELEMENT_RANGE = 1 << 16

def createIntMatrix(size:int) ->tuple[tuple[varDbl.VarDbl]]:
    '''
    Create an int matrix of size "size" with each element uniformly distributed between [-"ELEMENT_RANGE", +"ELEMENT_RANGE"]
    '''
    size = abs(size)
    return tuple([tuple([random.randint(-ELEMENT_RANGE, +ELEMENT_RANGE) for col in range(size)]) 
                  for row in range(size)])

def addNoise(ssMatrix:tuple[tuple[typing.Union[int, float]]], noise:float,
             retainIntValue:bool=True) -> tuple[tuple[varDbl.VarDbl]]:
    '''
    Add Gaussian noise of level "noise" to "ssMatrix", with the actual noise is adjusted by ELEMENT_RANGE.

    When "retainIntValue" is true, VarDbl.value() is int; Otherwise it is float.
    '''
    if not isSquareMatrix(ssMatrix, sType=(int, float)):
        raise ValueError(f'Invalid int or float')
    size = len(ssMatrix)
    noise *= ELEMENT_RANGE / math.sqrt(3)
    return tuple([tuple([varDbl.VarDbl(ssMatrix[row][col] + 
                                (int(noise * random.normalvariate()) if retainIntValue else noise * random.normalvariate()), noise) 
                    for col in range(size)]) for row in range(size)])


def determinant(ssMatrix:tuple[tuple[ElementType]]) -> ElementType:
    '''
    Calculate determinant and the adjugate matrix for "ssMatrix".
    If the ElementType contains  The result promotion is int -> Fraction -> float.
    '''
    isVar = isSquareMatrix(ssMatrix, sType=(varDbl.VarDbl,))
    if (not isSquareMatrix(ssMatrix, sType=(int, float, fractions.Fraction))) and (not isVar):
        raise ValueError(f'The input square matrix is illegal for determinant(): {ssMatrix}')
    size = len(ssMatrix)
    sPermut = {permut: permutSign(permut) for permut in itertools.permutations(range(size), size)}

    value = 0
    variance = 0
    for permut, sign in sPermut.items():
        val = sign
        if isVar:
            var = 1
            for x, y in enumerate(permut):
                val *= ssMatrix[x][y].value()
                var *= ssMatrix[x][y].variance()
            variance += var
        else:
            for x, y in enumerate(permut):
                val *= ssMatrix[x][y]
        value += val
    if not isVar:
        return value
    
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
                        var *= ssMatrix[x][y].variance()
                    else:
                        val *= ssMatrix[x][y].value()
                if sY in sVal:
                    sVal[sY] += val
                    assert sVar[sY] == var
                else:
                    sVal[sY] = val
                    sVar[sY] = var
        for sY, var in sVar.items():
            variance += (sVal[sY] ** 2) * var

    return varDbl.VarDbl(value, variance, True)



