'''
Extracting taylor expansion order as a numerical array {s1dTaylor} may ends up with 0 for higher-order expansions.

Thus, each function needs to be calculated independently

'''
import collections
import math
import logging
import typing
import unittest

import momentum
from varDbl import VarDbl, InitException

logger = logging.getLogger(__name__)


class Taylor1dException (Exception):
    def _name_() -> str:
        return 'Taylor1dException'

    def __init__(self, input:VarDbl, name:str, s1dTaylor:tuple[float], inPrec:bool, outPrec:bool, 
                 value:VarDbl, variance:VarDbl, n:int, newValue:VarDbl, newVariance:VarDbl, monotonics:int,
                 *args: object) -> None:
        super().__init__(*args)
        self.input = input
        self.name = name
        self.s1dTaylor = s1dTaylor
        self.inPrec = inPrec
        self.outPrec = outPrec

        self.value = value
        self.variance = variance
        self.n = n
        self.newValue = newValue
        self.newVariance = newVariance
        self.monotonics = monotonics

    def __str__(self) -> str:
        return (f'{self.__class__.__name__}: {self.name} for {self.input} at {self.n}:'
                f' val = {self.value}, var = {self.variance}, newVal = {self.newValue}, newVar = {self.newVariance},'
                f' monotonics={self.monotonics}')

class NotReliableException (Taylor1dException):
    pass

class NotPositiveException (Taylor1dException):
    pass

class NotMonotonicException (Taylor1dException):
    pass

class NotStableException (Taylor1dException):
    pass

class NotFiniteException (Taylor1dException):
    pass

Expansion = collections.namedtuple('Expansion', ('exp', 'val', 'var', 'limit', 'newVal', 'newVar', 'monotonics'))


class _Taylor:
    TAU = 7.18e-7   # The stability test 
    MIN_MONOTONIC_COUNT = 20

    __slots__ = ('_momentum')

    def __init__(self, binding:float = momentum.Momentum.BINDING_FACTOR) -> None:
        self._momentum = momentum.Momentum(binding=binding)


    @staticmethod
    def _writeList(fw, name, sList:tuple[float]):
        fw.write("Index:\t")
        for n in range(len(sList)):
            fw.write(f"{n}\t")
        fw.write("\n")
        fw.write(f"{name} Value:\t")
        hasVar = False
        for n in range(len(sList)):
            if isinstance(sList[n], VarDbl):    
                fw.write(f"{sList[n].value()}\t")
                hasVar = True
            else:
                fw.write(f"{sList[n]}\t")
        fw.write("\n")
        if hasVar:
            fw.write(f"{name} Uncertainty:\t")
            for n in range(len(sList)):
                if isinstance(sList[n], VarDbl):
                    fw.write(f"{sList[n].uncertainty()}\t")
                else:
                    fw.write("0\t")
            fw.write("\n")

    @staticmethod
    def headerForExpansion():
        return ("2n\tExponent Value\tExponent Variance\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty"
                "\tlimit\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\tMonotonic Count\n")

    @staticmethod
    def _writeExpansion(fw, n, varn, unc, value, variance, newValue, newVariance, monotonics):
        fw.write(f"{n}\t{varn.value()}\t{varn.variance()}"
                    f"\t{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}"
                    f"\t{unc}\t{newValue.value()}\t{newValue.variance()}\t{newVariance.value()}\t{newVariance.uncertainty()}"
                    f"\t{monotonics}\n")

    @staticmethod
    def headerForInput():
        return "name\tvalue\tuncertainty\tvariance\tinPrec\toutPrec\tBinding\tMaxOrder\n"
    
    @staticmethod
    def headerForResult():
        return "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\n"
    

    def _taylor1d(self, input:VarDbl, name:str, sTaylor:tuple[float], 
                 inPrec:bool, outPrec:bool, maxOrder:int=momentum.Momentum.MAX_ORDER, 
                 checkMonotonic=True, checkStability=True, checkReliablity=True, checkPositive=True,
                 dumpPath:str=None):
        '''
        1d Taylor expansion for differential series {s1dTaylor} at {input}, with {name} for logging.
        When {inPrec} is true, calculate Taylor expnasion against the precision of {input}.
        When {outPrec} is true, the result of the Taylor expnasion is the precision.
        s1dTaylor[n] should already normalized by /n!. 
        The max order of expansion is {maxOrder}, which should not exceed momentum.Momentum.MAX_ORDER

        When {checkMonotonic} is true, raise {NotMonotonicException} if 
            after full expansion, the monotonic count is still less than {MIN_MONOTONIC_COUNT}.
        It should always be True.

        When {checkStability} is true, raise {NotStableException} if 
            after full expansion, the value for the last expansion term is more than 
            TAU-fold of the expansion uncertainty.
        It should always be True.

        When {checkReliablity} is true, raise {NotReliableException} if
            the precision of the result variance is more than 1/5, in which 5 is the binding factor.
        It should always be True.

        When {checkPositive} is true, raise {NotPosive} if
            the expansion variance at any order becomes negative
        It should always be True.

        Both the result value and variance are guaranteed to be finite, otherwise
            raise {NotFiniteException} 

        Dump the expansion to {dumpPath} when it is provided.
        {dumpPath} can be read back and tested using verifyDumpFile()
        '''
        if (type(input) != VarDbl) or (not input.variance()):
            return VarDbl(sTaylor[0])
        n = 1
        for n in range(1, len(sTaylor)):
            if sTaylor[-n]:
                break
        s1dTaylor = sTaylor[:len(sTaylor) - n + 1]
        fw = None
        if dumpPath:
            fw = open(dumpPath, "w")
            fw.write(_Taylor.headerForInput())
            fw.write(f"{name}\t{input.value()}\t{input.uncertainty()}\t{input.variance()}"
                     f"\t{inPrec}\t{outPrec}"
                     f"\t{self._momentum._binding}\t{maxOrder}\n")
            _Taylor._writeList(fw, "Taylor", s1dTaylor)
            fw.write(_Taylor.headerForExpansion())

        monotonics = 0

        value = VarDbl(1, 0) if outPrec else s1dTaylor[0]
        variance = VarDbl()
        var = VarDbl(input.variance())
        if inPrec:
            var *= VarDbl( 1 /input.value() /input.value() )
        varn = VarDbl(var)
        prevVariance = VarDbl()
        for n in range(2, min(len(s1dTaylor)*2, maxOrder), 2):
            newValue = varn *s1dTaylor[n] *self._momentum[n] if n < len(s1dTaylor) else VarDbl()
            newVariance = 0
            for j in range(1, n):
                if (j < len(s1dTaylor)) and (n - j < len(s1dTaylor)):
                    try:
                        newVariance += varn * s1dTaylor[j] * s1dTaylor[n - j] * self._momentum[n]
                        if (n % 2) == 0:
                            newVariance -= varn * s1dTaylor[j] * self._momentum[j] * s1dTaylor[n - j] * self._momentum[n - j]
                    except (OverflowError, InitException) as ex:
                        if fw:
                            fw.write(f"NotFiniteException\t{ex}\t{n}\t{j}\t{s1dTaylor[j]}\t{s1dTaylor[n - j]}\n")
                            fw.close()
                        raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                                value, variance, n, newValue, newVariance, monotonics)

            if not isinstance(newVariance, VarDbl):
                newVariance = VarDbl(newVariance)
            try:
                value += newValue
                variance += newVariance
            except (OverflowError, InitException) as ex:
                if fw:
                    fw.write(f"NotFiniteException\t{ex}\t{n}\t{newVariance}\t{prevVariance}\t{variance}\n")
                    fw.close()
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
            unc = math.sqrt(variance.value()) *_Taylor.TAU if variance.value() > 0 else 0
            if abs(newVariance) <= abs(prevVariance):
                monotonics += 1
            else:
                monotonics = 0
            if fw:
                _Taylor._writeExpansion(fw, n, varn, unc, value, variance, newValue, newVariance, monotonics)

            if (not math.isfinite(value.value())) or (not math.isfinite(value.variance())) \
                    or (not math.isfinite(variance.value())) or (not math.isfinite(variance.variance())):
                if fw:
                    fw.write("NotFiniteException\n")
                    fw.close()
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
            if checkPositive and (variance.value() < 0):
                if fw:
                    fw.write("NotPositiveException\n")
                    fw.close()
                raise NotPositiveException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
            try:  
                varn *= var
            except InitException as ex:
                if fw:
                    fw.write(f"NotFiniteException\t{ex}\n")
                    fw.close()
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
            if varn.value() == 0:
                if fw:
                    fw.write(f"Expansion Truncation\n")
                break
            prevVariance = newVariance

        if checkMonotonic and (varn.value() > 0) and (monotonics < _Taylor.MIN_MONOTONIC_COUNT):
            if fw:
                fw.write("NotMonotonicException\n")
                fw.close()
            raise NotMonotonicException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        if checkStability and not ((abs(newValue.value()) < unc) or (abs(newValue.value()) < math.ulp(value.value()))):
            if fw:
                fw.write(f"NotStableException\t{unc}\n")
                fw.close()
            raise NotStableException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        
        if outPrec:
            try:
                value *= s1dTaylor[0]
                variance *= s1dTaylor[0]**2
            except InitException as ex:
                if fw:
                    fw.write(f"NotFiniteException\t{ex}\n")
                    fw.close()
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
        if checkReliablity and (variance.value() * self._momentum.binding < variance.uncertainty()):
            if fw:
                fw.write("NotReliableException\n")
                fw.close()
            raise NotReliableException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        if fw:
            fw.write(_Taylor.headerForResult())
            fw.write(f"{value.value()}\t{value.variance()}\t{variance.value()}\t{variance.variance()}\n")
            fw.close()
        try:
            if variance:
                return VarDbl(value.value(), variance.value() + value.variance(), True) 
            else:
                return value
        except InitException:
            if fw:
                fw.write("NotFiniteException\n")
                fw.close()
            raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance)

    @staticmethod   
    def verifyDumpFile(self:unittest.TestCase, dumpPath:str, checkExp=True)\
            -> tuple[VarDbl, list[float], Expansion, typing.Union[VarDbl, str]]:
        '''
        When input uncertainty is rounding error, checkExp=False
        '''
        with open(dumpPath) as f:
            self.assertEqual(next(f), _Taylor.headerForInput())
            sWord = next(f).strip().split('\t')
            x = VarDbl(sWord[1], sWord[2])
            self.assertEqual(len(sWord), len(_Taylor.headerForInput().split('\t')))
            sWord = next(f).strip().split('\t')
            self.assertEqual(sWord[0], 'Index:')
            n = len(sWord) - 1
            self.assertListEqual(sWord[1:], [f'{i}' for i in range(n)])
            sWord = next(f).strip().split('\t')
            self.assertEqual(len(sWord), n + 1)
            self.assertEqual(sWord[0], 'Taylor Value:')
            sTaylor = list(map(float, sWord[1:]))
            self.assertEqual(next(f), _Taylor.headerForExpansion())
            sExpansion = []
            infinite = False
            for n in range(2, momentum.Momentum.MAX_ORDER, 2):
                sWord = next(f).strip().split('\t')
                try:
                    nn = int(sWord[0])
                except ValueError:
                    break
                self.assertEqual(n, nn)
                try:
                    exp = VarDbl(sWord[1], sWord[2])
                    if checkExp:
                        self.assertAlmostEqual(x.uncertainty() **n / exp.value(), 1)
                    val = VarDbl(sWord[3], sWord[4])
                    var = VarDbl(sWord[5], sWord[6])
                    limit = float(sWord[7])
                    if limit:
                        self.assertAlmostEqual(math.sqrt(var.value()) * _Taylor.TAU / limit, 1)
                    else:
                        self.assertLessEqual(var.value(), 0)
                    newVal = VarDbl(sWord[8], sWord[9])
                    newVar = VarDbl(sWord[10], sWord[11])
                except InitException:
                    infinite = True
                monotonics = int(sWord[12])
                sExpansion.append(Expansion(exp, val, var, limit, newVal, newVar, monotonics))
            else:
                sWord = next(f).strip().split('\t')
            if _Taylor.headerForResult().startswith(sWord[0]):
                sWord = next(f).strip().split('\t')
                res = VarDbl(sWord[0], float(sWord[1]) + float(sWord[2]), True)
            else:
                res = sWord[0]
                if infinite:
                    self.assertEqual(res, 'NotFiniteException')
            return x, sTaylor, sExpansion, res
    


class Taylor:  
    TAU = _Taylor.TAU
    MIN_MONOTONIC_COUNT = _Taylor.MIN_MONOTONIC_COUNT
    MAX_POLY_ORDER = momentum.Momentum.MAX_ORDER // 2
    _taylor = _Taylor()

    @staticmethod
    def maxOrder():
        return momentum.Momentum.MAX_ORDER
    
    @staticmethod
    def taylor1d(input:VarDbl, name:str, sTaylor:tuple[float], 
                 inPrec:bool, outPrec:bool, maxOrder:int=momentum.Momentum.MAX_ORDER,
                 dumpPath:str=None):
        return Taylor._taylor._taylor1d(input, name, sTaylor, 
                 inPrec, outPrec, maxOrder=maxOrder,
                 dumpPath=dumpPath)
    
    @staticmethod   
    def verifyDumpFile(self:unittest.TestCase, dumpPath:str, checkExp=True):
        return _Taylor.verifyDumpFile(self, dumpPath, checkExp=checkExp)

    @staticmethod
    def polynominal(input:VarDbl, sCoeff:tuple[float],
                    dumpPath:str=None):
        '''
        1d Taylor expansion for polynominal at "input" with "sCoeff".
        Allow input.value() +- input.uncertainty() to include 0
        '''
        if len(sCoeff) > (Taylor.MAX_POLY_ORDER):
            raise ValueError(f'The lenght {len(sCoeff)} of polynominal coefficient is more than half of {Taylor.maxOrder()}: {sCoeff}')
        exp = len(sCoeff) - 1
        s1dTaylor = [sCoeff[0]] + [0] * 2*exp
        sPow = [1, input.value()]

        for j, coeff in enumerate(sCoeff):
            if not j:
                continue
            sTaylor = [1, j]
            for k in range(2, j + 1):
                sTaylor.append( sTaylor[-1] * (j + 1 - k)//k )
                while k >= len(sPow):
                    sPow.append(sPow[k-1] * input.value())
            for k in range(j + 1):
                s1dTaylor[k] += coeff * sTaylor[k] * sPow[j - k]
                 
        return Taylor._taylor._taylor1d(input, f'poly({sCoeff})', s1dTaylor, False, False, 
                dumpPath = dumpPath, checkMonotonic = False, checkStability = False)

    @staticmethod
    def exp(input:VarDbl, dumpPath:str=None) -> VarDbl:
        sTaylor = [math.exp(input.value()), 1.0]
        for i in range(2, Taylor.maxOrder()):
            sTaylor.append(sTaylor[-1]/i)
            if not sTaylor[-1]:
                break
        return Taylor._taylor._taylor1d(input, f"exp({input})", sTaylor, False, True, 
                    dumpPath=dumpPath)
    
    @staticmethod
    def log(input:VarDbl, dumpPath:str=None) -> VarDbl:
        sTaylor = []
        sTaylor.append(math.log(input.value()))
        for i in range(1, Taylor.maxOrder()):
            sTaylor.append(1/i if ((i%2) == 1) else -1/i)
        return Taylor._taylor._taylor1d(input, f"log({input})", sTaylor, True, False, 
                    dumpPath=dumpPath)

    @staticmethod
    def sin(input:VarDbl, dumpPath:str=None) -> VarDbl:
        sTaylor = []
        x = input.value()
        sTaylor.append( math.sin(x) )
        n = 1.0
        for i in range(1, Taylor.maxOrder()):
            n *= i
            if not math.isfinite(n):
                break
            match i % 4:
                case 0:
                    sTaylor.append(math.sin(x) / n)
                case 1:
                    sTaylor.append(math.cos(x) / n)  
                case 2:
                    sTaylor.append(-math.sin(x) / n)
                case 3:
                    sTaylor.append(-math.cos(x) / n)
        return Taylor._taylor._taylor1d(input, f"sin({input})", sTaylor, False, False, 
                    dumpPath=dumpPath)
    
    @staticmethod
    def pow(input:VarDbl, exp:float, dumpPath:str=None) -> VarDbl:
        match exp:
            case 0:
                return VarDbl(1, 0)
            case 1:
                return VarDbl(input)
        if (exp > 0) and ((type(exp) == int) or (math.ceil(exp) == math.floor(exp))):
            sCoeff = [0] * int(exp)
            sCoeff.append(1)
            return Taylor.polynominal(input, sCoeff, dumpPath=dumpPath)
        sTaylor = [math.pow(input.value(), exp), exp]
        for i in range(2, Taylor.maxOrder()):
            sTaylor.append( sTaylor[-1]/i * (exp + 1 - i) )
            if not sTaylor[-1]:
                break
        return Taylor._taylor._taylor1d(input, f"({input})**{exp}", sTaylor, True, True, 
                    dumpPath=dumpPath)
    

