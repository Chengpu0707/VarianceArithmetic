'''
Extracting taylor expansion order as a numerical array {s1dTaylor} may ends up with 0 for higher-order expansions.

Thus, each function needs to be calculated independently

'''
import collections
import math
import logging
import scipy.special
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

Expansion = collections.namedtuple('Expansion', ('exp', 'val', 'var', 'newVal', 'newVar', 'monotonics'))


class Taylor:
    MIN_MONOTONIC_COUNT = 20

    @staticmethod
    def default():
        return _taylor

    __slots__ = ('_momentum', '_leakage')

    def __init__(self, bounding:float = 5.0) -> None:
        self._momentum = momentum.Normal(bounding=bounding)
        self._leakage = 1 - scipy.special.erf(bounding/math.sqrt(2))

    @property
    def momentum(self):
        return self._momentum

    @property
    def leakage(self):
        return self._leakage   

    @staticmethod
    def _writeList(fw, name, sList:tuple[typing.Union[float, VarDbl]]):
        fw.write("Index:\t")
        for n in range(len(sList)):
            fw.write(f"{n}\t")
        fw.write("\n")
        fw.write(f"{name} Value:\t")
        hasVar = False
        for n in range(len(sList)):
            if isinstance(sList[n], VarDbl):    
                fw.write(f"{sList[n].value()}\t")
                if sList[n].variance() > 0:
                    hasVar = True
            else:
                fw.write(f"{sList[n]}\t")
        fw.write("\n")
        if hasVar:
            fw.write(f"{name} Uncertainty:\t")
            for n in range(len(sList)):
                fw.write(f"{sList[n].uncertainty()}\t")
            fw.write("\n")

    @staticmethod
    def headerForExpansion():
        return ("2n\tMonotonics\tExponent\tMomentum\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty"  # 0-5
                "\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\n")

    @staticmethod
    def headerForInput():
        return ("name\tvalue\tuncertainty\tvariance\tinPrec\toutPrec\tBinding\tMaxOrder\tMinMonotonic"
                "\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive\tenableExpansionTruncation\n")
    
    @staticmethod
    def headerForResult():
        return "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\n"
    

    def taylor1d(self, input:VarDbl, name:str, s1dTaylor:tuple[typing.Union[float, VarDbl]], 
                 inPrec:bool, outPrec:bool, maxOrder:int=0, 
                 checkMinMonotonic=True, checkStability=True, checkReliablity=True, checkPositive=True,
                 dumpPath:str=None):
        '''
        1d Taylor expansion for differential series {s1dTaylor} at {input}, with {name} for logging.
        When {inPrec} is true, calculate Taylor expnasion against the precision of {input}.
        When {outPrec} is true, the result of the Taylor expnasion is the precision.
        s1dTaylor[n] should already normalized by /n!. 
        The max order of expansion is {self.momentum.maxOrder}, which should not exceed {self.momentum.maxOrder}

        When {checkMinMonotonic} is true, raise {NotMonotonicException} if 
            after full expansion, the monotonic count is still less than {MIN_MONOTONIC_COUNT}.
        When the momentum is more than {MAX_MONOTONIC_MOMENTUM}, turn off resetting the monotonic check
            due to numerical unstability, e.g., TestDumpFile.test_NotMonotonicException_AfterMax().
        It should always be True.

        When {checkStability} is true, raise {NotStableException} if 
            after full expansion, the value for the last expansion term is more than 
            TAU-fold of the expansion uncertainty.
        It should always be True.

        When {checkReliablity} is true, raise {NotReliableException} if
            the precision of the result variance is more than the bounding factor.
        It should always be True.

        When {checkPositive} is true, raise {NotPosive} if
            the expansion variance at any order becomes negative
        It should always be True.

        Both the result value and variance are guaranteed to be finite, otherwise
            raise {NotFiniteException}
        But if the monotonic count is more than {MIN_MONOTONIC_COUNT},
            roll back to the previous value and variance.

        Dump the expansion to {dumpPath} when it is provided.
        {dumpPath} can be read back and tested using verifyDumpFile()
        '''
        if not maxOrder:
            maxOrder = self._momentum.maxOrder
        for n in range(len(s1dTaylor)):
            if isinstance(s1dTaylor[n], VarDbl):
                if (not math.isfinite(s1dTaylor[n].value())) or (not math.isfinite(s1dTaylor[n].variance())):
                    raise ValueError(f'Taylor [{n}]={s1dTaylor[n]}')
            elif isinstance(s1dTaylor[n], float):
                if not math.isfinite(s1dTaylor[n]):
                    raise ValueError(f'Taylor [{n}]={s1dTaylor[n]}')
            elif isinstance(s1dTaylor[n], int):
                pass
            else:
                raise ValueError(f'Taylor [{n}]={s1dTaylor[n]}')
            
        if (type(input) != VarDbl) or (not input.variance()):
            return VarDbl(s1dTaylor[0])
        fw = None
        if dumpPath:
            fw = open(dumpPath, "w")
            fw.write(Taylor.headerForInput())
            fw.write(f"{name}\t{input.value()}\t{input.uncertainty()}\t{input.variance()}"
                     f"\t{inPrec}\t{outPrec}"
                     f"\t{self._momentum.bounding}\t{maxOrder}\t{Taylor.MIN_MONOTONIC_COUNT}"
                     f"\t{checkMinMonotonic}\t{checkStability}\t{checkReliablity}\t{checkPositive}\tFalse\n")
            Taylor._writeList(fw, "Taylor1d", s1dTaylor)
            fw.write(Taylor.headerForExpansion())

        monotonics = 0
        monotonicPrev = True

        value = VarDbl(1, 0) if outPrec else VarDbl(s1dTaylor[0])
        variance = VarDbl()
        var = input.variance()
        if inPrec:
            var *= 1 /input.value() /input.value()
        varn = var
        prevValue = VarDbl()
        prevVariance = VarDbl()
        for n in range(2, min(len(s1dTaylor), maxOrder), 2):
            oldValue = VarDbl(value)
            oldVariance = VarDbl(variance)
            finite = True
            try:
                newValue = s1dTaylor[n] * varn * self._momentum[n]
                newVariance = 0
                for j in range(1, n):
                    newVariance += s1dTaylor[j] * s1dTaylor[n - j] * \
                                    varn * (self._momentum[n] - self._momentum[j] * self._momentum[n - j])
                if not isinstance(newValue, VarDbl):
                    newValue = VarDbl(newValue)
                if not isinstance(newVariance, VarDbl):
                    newVariance = VarDbl(newVariance)
                value += newValue
                variance += newVariance
                if (not math.isfinite(value.value())) or (not math.isfinite(value.variance() + variance.value())) or \
                        (not math.isfinite(variance.variance())):
                    finite = False
                elif abs(newVariance.value()) <= abs(prevVariance.value()):
                    monotonics += 1
                elif (monotonics >= Taylor.MIN_MONOTONIC_COUNT) and monotonicPrev:
                    monotonicPrev = False
                else:
                    monotonics = 0
                prevValue = newValue
                prevVariance = newVariance

                if fw:
                    fw.write(f"{n}\t{monotonics}\t{varn}\t{self._momentum[n]}"
                             f"\t{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}"
                             f"\t{newValue.value()}\t{newValue.variance()}\t{newVariance.value()}\t{newVariance.uncertainty()}"
                             f"\n")
            except (OverflowError, InitException) as ex:
                finite = False
            except BaseException as ex:
                if fw:
                    fw.write(f"Fail\t{ex}\t{n}\t{j}\t{s1dTaylor[j]}\t{s1dTaylor[n - j]}\t{self._momentum[n]}\n")
                raise ex
            if not finite:
                if checkMinMonotonic and (monotonics >= Taylor.MIN_MONOTONIC_COUNT):
                    # roll back to the previous value
                    value = oldValue
                    variance = oldVariance
                    newValue = prevValue
                    newVariance = prevVariance
                    break
                if not math.isfinite(value.value()):
                    if fw:
                        fw.write("NotFiniteException\tvalue\n")
                        fw.close()
                elif (not math.isfinite(value.variance() + variance.value())) or (not math.isfinite(variance.variance())):
                    if fw:
                        fw.write("NotFiniteException\tvariance\n")
                        fw.close()
                else:
                    if fw:
                        fw.write(f"NotFiniteException\tcalculation\t{n}\t{j}\t{value}\t{newValue}\t{variance}\t{newVariance}\n")
                        fw.close()
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
                
            if checkPositive and (variance.value() < 0):
                if fw:
                    fw.write("NotPositiveException\n")
                    fw.close()
                raise NotPositiveException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
            varn *= var
            if not math.isfinite(varn):
                break
            if varn == 0:
                break

        if checkMinMonotonic and (varn > 0) and (monotonics < Taylor.MIN_MONOTONIC_COUNT):
            if fw:
                fw.write("NotMonotonicException\n")
                fw.close()
            raise NotMonotonicException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        unc = math.sqrt(value.variance() + variance.value()) * self.leakage
        if checkStability and not ((abs(newValue.value()) < unc) or (abs(newValue.value()) < math.ulp(value.value()))):
            if fw:
                fw.write(f"NotStableException\t{n}\t{len(s1dTaylor)}\t{unc}\n")
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
        if checkReliablity and (variance.value() * self._momentum.bounding < variance.uncertainty()):
            if fw:
                fw.write("NotReliableException\n")
                fw.close()
            raise NotReliableException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        if fw:
            fw.write(Taylor.headerForResult())
            fw.write(f"{value.value()}\t{value.variance()}\t{variance.value()}\t{variance.variance()}\n")
            fw.close()
        try:
            return VarDbl(value.value(), variance.value() + value.variance(), True) 
        except InitException as ex:
            if fw:
                fw.write(f"NotFiniteException\tresult\t{ex}\n")
                fw.close()
            raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance)


    @staticmethod   
    def verifyDumpFile(testcase:unittest.TestCase, dumpPath:str, checkExp=True)\
            -> tuple[VarDbl, list[float], Expansion, typing.Union[VarDbl, str]]:
        '''
        When input uncertainty is rounding error, checkExp=False
        '''
        with open(dumpPath) as f:
            testcase.assertEqual(next(f), Taylor.headerForInput())
            sWord = next(f).strip().split('\t')
            x = VarDbl(sWord[1], sWord[2])
            testcase.assertEqual(len(sWord), len(Taylor.headerForInput().split('\t')))
            sWord = next(f).strip().split('\t')
            testcase.assertEqual(sWord[0], 'Index:')
            n = len(sWord) - 1
            testcase.assertListEqual(sWord[1:], [f'{i}' for i in range(n)])
            sWord = next(f).strip().split('\t')
            testcase.assertEqual(len(sWord), n + 1)
            testcase.assertEqual(sWord[0], 'Taylor1d Value:')
            sTaylor = list(map(float, sWord[1:]))
            testcase.assertEqual(next(f), Taylor.headerForExpansion())
            sExpansion = []
            infinite = False
            mmt = momentum.Normal()
            for n in range(2, mmt.maxOrder, 2):
                sWord = next(f).strip().split('\t')
                try:
                    nn = int(sWord[0])
                except ValueError:
                    break
                testcase.assertEqual(n, nn)
                try:
                    monotonics = int(sWord[1])
                    exp = float(sWord[2])
                    if checkExp:
                        testcase.assertAlmostEqual(x.uncertainty() **n / exp, 1)
                    val = VarDbl(sWord[4], sWord[5])
                    var = VarDbl(sWord[6], sWord[7])
                    newVal = VarDbl(sWord[8], sWord[9])
                    newVar = VarDbl(sWord[10], sWord[11])
                except InitException:
                    infinite = True
                sExpansion.append(Expansion(exp, val, var, newVal, newVar, monotonics))
            else:
                sWord = next(f).strip().split('\t')
            if Taylor.headerForResult().startswith(sWord[0]):
                sWord = next(f).strip().split('\t')
                res = VarDbl(sWord[0], float(sWord[1]) + float(sWord[2]), True)
            else:
                res = sWord[0]
                if infinite:
                    testcase.assertEqual(res, 'NotFiniteException')
            return x, sTaylor, sExpansion, res
    


    def polynominal(self, input:VarDbl, sCoeff:tuple[float],
                    dumpPath:str=None):
        '''
        1d Taylor expansion for polynominal at "input" with "sCoeff".
        Allow input.value() +- input.uncertainty() to include 0
        '''
        if len(sCoeff) > (self.momentum.maxOrder // 2):
            raise ValueError(f'The lenght {len(sCoeff)} of polynominal coefficient is more than half of {self.momentum.maxOrder // 2}: {sCoeff}')
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
                 
        return self.taylor1d(input, f'poly({sCoeff})', s1dTaylor, False, False, 
                dumpPath = dumpPath, checkMinMonotonic = False, checkStability = False)

    def exp(self, input:VarDbl, dumpPath:str=None) -> VarDbl:
        sTaylor = [math.exp(input.value()), 1.0]
        for i in range(2, self.momentum.maxOrder):
            sTaylor.append(sTaylor[-1]/i)
        return self.taylor1d(input, f"exp({input})", sTaylor, False, True, 
                    dumpPath=dumpPath)
    
    def log(self, input:VarDbl, dumpPath:str=None) -> VarDbl:
        sTaylor = []
        sTaylor.append(math.log(input.value()))
        for i in range(1, self.momentum.maxOrder):
            sTaylor.append(1/i if ((i%2) == 1) else -1/i)
        return self.taylor1d(input, f"log({input})", sTaylor, True, False, 
                    dumpPath=dumpPath)

    def sin(self, input:VarDbl, dumpPath:str=None) -> VarDbl:
        sTaylor = []
        x = input.value()
        sTaylor.append( math.sin(x) )
        fac = 1.0
        for i in range(1, self.momentum.maxOrder):
            fac /= i
            match i % 4:
                case 0:
                    sTaylor.append(math.sin(x) *fac)
                case 1:
                    sTaylor.append(math.cos(x) *fac)  
                case 2:
                    sTaylor.append(-math.sin(x) *fac)
                case 3:
                    sTaylor.append(-math.cos(x) *fac)
        return self.taylor1d(input, f"sin({input})", sTaylor, False, False, 
                    dumpPath=dumpPath)
    
    def pow(self, input:VarDbl, exp:float, dumpPath:str=None) -> VarDbl:
        match exp:
            case 0:
                return VarDbl(1, 0)
            case 1:
                return VarDbl(input)
        if (exp > 0) and ((type(exp) == int) or (math.ceil(exp) == math.floor(exp))):
            sCoeff = [0] * int(exp)
            sCoeff.append(1)
            return self.polynominal(input, sCoeff, dumpPath=dumpPath)
        sTaylor = [math.pow(input.value(), exp), exp]
        for i in range(2, self.momentum.maxOrder):
            sTaylor.append( sTaylor[-1]/i * (exp + 1 - i) )
        return self.taylor1d(input, f"({input})**{exp}", sTaylor, True, True, 
                    dumpPath=dumpPath)
    

_taylor = Taylor()  # Singleton instance   
