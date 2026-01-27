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

Expansion = collections.namedtuple('Expansion', ('order', 'taylor', 'exp', 'momentum', 'monotonics', 'val', 'var', 'newVal', 'newVar'))



class Taylor:
    def __new__(cls):
        raise TypeError('Static classes cannot be instantiated')
    
    MIN_MONOTONIC_COUNT = 20

    DUMP_PATH_INPUT_HEADER = (
        "result\tvalue\tuncertainty\tinPrec\toutPrec\tbounding\tmaxOrder\tMinMonotonic"
        "\tcheckMonotonic\tcheckStability\tcheckReliablity\tcheckPositive\tName\n")
    DUMP_PATH_EXPANSION_HEADER = (
        "Order\tTaylor Value\tTaylor Uncertainty\tExponent\tMomentum\tMonotonics"
        "\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty"
        "\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\n")
    DUMP_PATH_OUTPUT_HEADER = (
        "Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\tException\n")  

    @staticmethod
    def _writeList(fw, name, sList:tuple[typing.Union[float, VarDbl]]):
        fw.write(f"{name}\t")
        for n in range(len(sList)):
            fw.write(f"{n}\t")
        fw.write("\n")
        hasVar = False
        for n in range(len(sList)):
            if isinstance(sList[n], VarDbl):    
                if sList[n].variance() > 0:
                    hasVar = True
        fw.write("Value\t" if hasVar else "Double\t")
        for n in range(len(sList)):
            if isinstance(sList[n], VarDbl):    
                fw.write(f"{sList[n].value()}\t")
            else:
                fw.write(f"{sList[n]}\t")
        fw.write("\n")
        if hasVar:
            fw.write("Uncertainty\t")
            for n in range(len(sList)):
                fw.write(f"{sList[n].uncertainty()}\t")
            fw.write("\n")

    @staticmethod
    def _writeResult(fw, value:VarDbl, variance:VarDbl, exception:str):
        if not fw:
            return
        fw.write(Taylor.DUMP_PATH_OUTPUT_HEADER)
        fw.write(f"{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}\t{exception}\n")
        fw.close()
    

    @staticmethod
    def taylor1d(input:VarDbl, name:str, s1dTaylor:tuple[typing.Union[float, VarDbl]], 
                 inPrec:bool, outPrec:bool, 
                 momentum=momentum.NORMAL,
                 checkMinMonotonic=True, checkStability=True, checkReliablity=True, checkPositive=True, checkLSB=False,
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
            fw.write(Taylor.DUMP_PATH_INPUT_HEADER)
            fw.write(f"{s1dTaylor[0]}\t{input.value()}\t{input.uncertainty()}\t{inPrec}\t{outPrec}"
                     f"\t{momentum.bounding}\t{momentum.maxOrder}\t{Taylor.MIN_MONOTONIC_COUNT}"
                     f"\t{checkMinMonotonic}\t{checkStability}\t{checkReliablity}\t{checkPositive}"
                     f"\t{name}\n")
            fw.write(Taylor.DUMP_PATH_EXPANSION_HEADER)

        monotonics = 0
        monotonicPrev = True

        value = VarDbl(1, 0) if outPrec else VarDbl(s1dTaylor[0])
        variance = VarDbl()
        unc = input.uncertainty()
        if inPrec:
            unc /= input.value()
        uncN = 1
        prevVariance = VarDbl()
        infinite = None
        for n in range(1, min(len(s1dTaylor), momentum.maxOrder)):
            uncN *= unc
            if not math.isfinite(uncN):
                break
            if uncN == 0:
                break
            try:
                if abs(uncN) < 1:
                    newValue = s1dTaylor[n] * (uncN * momentum[n])
                else:
                    newValue = s1dTaylor[n] * uncN * momentum[n]
                if not isinstance(newValue, VarDbl):
                    newValue = VarDbl(newValue)

                try:
                    newVariance = 0
                    for j in range(1, n):
                        if abs(uncN) < 1:
                            newVariance += s1dTaylor[j] * s1dTaylor[n - j] * (uncN * \
                                            (momentum[n] - momentum[j] * momentum[n - j]))
                        else:
                            newVariance += s1dTaylor[j] * s1dTaylor[n - j] * uncN * \
                                            (momentum[n] - momentum[j] * momentum[n - j])
                    if not isinstance(newVariance, VarDbl):
                        newVariance = VarDbl(newVariance)
                except (OverflowError, InitException) as ex:
                    infinite = 'NotFiniteException\tnewVariance'

                try:    
                    value += newValue
                except (OverflowError, InitException) as ex:
                    infinite = 'NotFiniteException\tvalue'
                try:
                    variance += newVariance
                except (OverflowError, InitException) as ex:
                    infinite = 'NotFiniteException\tvariance'

                if (not infinite) and ((n & 1) == 0):
                    if abs(newVariance.value()) <= abs(prevVariance.value()):
                        monotonics += 1
                    elif (monotonics >= Taylor.MIN_MONOTONIC_COUNT) and monotonicPrev:
                        monotonicPrev = False
                    else:
                        monotonics = 0
                    prevVariance = newVariance

                if fw:
                    tyVal = s1dTaylor[n].value() if type(s1dTaylor[n]) == VarDbl else s1dTaylor[n]
                    tyVar = s1dTaylor[n].uncertainty() if type(s1dTaylor[n]) == VarDbl else 0
                    fw.write(f"{n}\t{tyVal}\t{tyVar}\t{uncN}\t{momentum[n]}\t{monotonics}"
                             f"\t{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}"
                             f"\t{newValue.value()}\t{newValue.variance()}\t{newVariance.value()}\t{newVariance.uncertainty()}"
                             f"\n")
            except (OverflowError, InitException) as ex:
                infinite = 'NotFiniteException\tnewValue'
            except BaseException as ex:
                infinite = f'NotFiniteException\t{ex}'
                if fw:
                    fw.write(f"{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}\t{ex}\n")
                raise ex
            if infinite:
                Taylor._writeResult(fw, value, variance, infinite)
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
                
            if checkPositive and (variance.value() < 0):
                Taylor._writeResult(fw, value, variance, 'NotPositiveException')
                raise NotPositiveException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)

        if checkMinMonotonic and (uncN > 0) and (monotonics < Taylor.MIN_MONOTONIC_COUNT):
            Taylor._writeResult(fw, value, variance, 'NotMonotonicException')
            raise NotMonotonicException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        if checkLSB:
            if abs(newValue.value()) >= math.ulp(value.value()):
                Taylor._writeResult(fw, value, variance, 'NotStableException\tvalue')
                raise NotStableException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
        if checkStability:
            unc = math.sqrt(value.variance() + variance.value()) * momentum.leakage
            if unc > 0 and abs(newValue.value()) >= unc:
                Taylor._writeResult(fw, value, variance, 'NotStableException\tuncertainty')
                raise NotStableException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
        
        if outPrec:
            try:
                value *= s1dTaylor[0]
            except InitException:
                Taylor._writeResult(fw, value, variance, f'NotFiniteException\toutPrec value\t{s1dTaylor[0]}')
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
            try:
                variance *= s1dTaylor[0]**2
            except InitException:
                Taylor._writeResult(fw, value, variance, f'NotFiniteException\toutPrec variance\t{s1dTaylor[0]}')
                raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance, monotonics)
        if checkReliablity and (variance.value() * momentum.bounding < variance.uncertainty()):
            Taylor._writeResult(fw, value, variance, f'NotReliableException')
            raise NotReliableException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance, monotonics)
        try:
            ret = VarDbl(value.value(), math.sqrt(variance.value() + value.variance()))
        except InitException as ex:
            Taylor._writeResult(fw, value, variance, f'NotFiniteException\t{ex}')
            raise NotFiniteException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance)
        Taylor._writeResult(fw, value, variance, '')
        return ret


    @staticmethod   
    def verifyDumpFile(testcase:unittest.TestCase, dumpPath:str)\
            -> tuple[VarDbl, list[Expansion], typing.Union[VarDbl, str]]:
        '''
        When input uncertainty is rounding error, checkExp=False
        '''
        with open(dumpPath) as f:
            testcase.assertEqual(next(f), Taylor.DUMP_PATH_INPUT_HEADER)
            sHdr = Taylor.DUMP_PATH_INPUT_HEADER.split('\t')
            sWord = next(f).strip().split('\t')
            testcase.assertEqual(len(sWord), len(sHdr))
            for i in (8,9,10,11):
                testcase.assertIn(sWord[i], ('True', 'true', '1'))
            for i in (3,4):
                testcase.assertIn(sWord[i], ('True', 'true', '1', 'False', 'false', '0'))
            testcase.assertEqual(int(sWord[7]), Taylor.MIN_MONOTONIC_COUNT)
            sInput = {
                'result': float(sWord[0]),
                'input': VarDbl(sWord[1], sWord[2]),
                'inPrec': sWord[3] in ('True', 'true', '1'), 'outPrec': sWord[4] in ('True', 'true', '1'),
                'bounding': float(sWord[5]), 'maxOrder': int(sWord[6]),
             }

            testcase.assertEqual(next(f), Taylor.DUMP_PATH_EXPANSION_HEADER)
            sExpansion = [0]
            infinite = False
            size = len(Taylor.DUMP_PATH_EXPANSION_HEADER.split('\t'))
            for n in range(1, 10000):
                try:
                    sWord = next(f).strip().split('\t')
                except StopIteration:
                    testcase.fail(f'Unexpected end of file {dumpPath} at line {n+1}: {sWord}')
                    break
                if size != len(sWord):
                    break
                nn = int(sWord[0])
                testcase.assertEqual(n, nn)
                try:
                    ty = VarDbl(sWord[1], sWord[2])
                    exp, mmt = map(float, sWord[3:5])
                    monotonics = int(sWord[5])
                    if n & 1:
                        testcase.assertEqual(0, mmt)
                    else:
                        testcase.assertLess(0, mmt)
                    val = VarDbl(sWord[6], sWord[7])
                    var = VarDbl(sWord[8], sWord[9])
                    newVal = VarDbl(sWord[10], sWord[11])
                    newVar = VarDbl(sWord[12], sWord[13])
                except InitException:
                    infinite = True
                sExpansion.append(Expansion(n, ty, exp, mmt, monotonics, val, var, newVal, newVar))
            testcase.assertListEqual(sWord, Taylor.DUMP_PATH_OUTPUT_HEADER.strip().split('\t'))
            sWord = next(f).strip().split('\t')
            testcase.assertLessEqual(4, len(sWord))
            if infinite:
                testcase.assertEqual(sWord[4], 'NotFiniteException')
            if 4 < len(sWord):
                res = sWord[4]
            else:
                val = VarDbl(sWord[0], sWord[1])
                var = VarDbl(sWord[2], sWord[3])
                res = VarDbl(val.value(), math.sqrt(val.variance() + var.value()))
            return sInput, sExpansion, res
    


    @staticmethod
    def polynominal1d(input:VarDbl, sCoeff:tuple[float], 
                      momentum=momentum.NORMAL,
                      dumpPath:str=None):
        '''
        1d Taylor expansion for polynominal at "input" with "sCoeff".
        Allow input.value() +- input.uncertainty() to include 0
        '''
        if len(sCoeff) > (momentum.maxOrder // 2):
            raise ValueError(f'The lenght {len(sCoeff)} of polynominal coefficient is more than half of {momentum.maxOrder}: {sCoeff}')
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
                 
        return Taylor.taylor1d(input, f'poly({sCoeff})', s1dTaylor, False, False, 
                momentum = momentum,
                dumpPath = dumpPath, checkMinMonotonic = False, checkStability = False)

    @staticmethod
    def exp(input:VarDbl, momentum=momentum.NORMAL, dumpPath:str=None) -> VarDbl:
        sTaylor = [math.exp(input.value()), 1.0]
        for i in range(2, momentum.maxOrder):
            sTaylor.append(sTaylor[-1]/i)
        return Taylor.taylor1d(input, f"exp({input})", sTaylor, False, True, 
                               momentum=momentum, dumpPath=dumpPath)
    
    @staticmethod
    def log(input:VarDbl, momentum=momentum.NORMAL, dumpPath:str=None) -> VarDbl:
        sTaylor = []
        sTaylor.append(math.log(input.value()))
        for i in range(1, momentum.maxOrder):
            sTaylor.append(1/i if ((i%2) == 1) else -1/i)
        return Taylor.taylor1d(input, f"log({input})", sTaylor, True, False, 
                               momentum=momentum, dumpPath=dumpPath)

    @staticmethod
    def sin(input:VarDbl, momentum=momentum.NORMAL, dumpPath:str=None) -> VarDbl:
        sTaylor = []
        x = input.value()
        sTaylor.append( math.sin(x) )
        fac = 1.0
        for i in range(1, momentum.maxOrder):
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
        return Taylor.taylor1d(input, f"sin({input})", sTaylor, False, False, 
                               momentum=momentum, dumpPath=dumpPath)
    
    @staticmethod
    def pow(input:VarDbl, exp:float, momentum=momentum.NORMAL, dumpPath:str=None) -> VarDbl:
        match exp:
            case 0:
                return VarDbl(1, 0)
            case 1:
                return VarDbl(input)
        if (exp > 0) and ((type(exp) == int) or (math.ceil(exp) == math.floor(exp))):
            sCoeff = [0] * int(exp)
            sCoeff.append(1)
            return Taylor.polynominal1d(input, sCoeff, momentum=momentum, dumpPath=dumpPath)
        sTaylor = [math.pow(input.value(), exp), VarDbl(exp)]
        for i in range(2, momentum.maxOrder):
            sTaylor.append( sTaylor[-1] * ((exp + 1 - i)/i) )
        return Taylor.taylor1d(input, f"({input})**{exp}", sTaylor, True, True, 
                               momentum=momentum, dumpPath=dumpPath)
    

