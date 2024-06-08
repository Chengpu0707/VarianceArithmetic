
import math
import logging
import typing

import varDbl
import momentum

logger = logging.getLogger(__name__)


class NotReliableException (Exception):
    def __init__(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[varDbl.VarDbl], inPrec:bool, outPrec:bool,
                 value:varDbl.VarDbl, variance:varDbl.VarDbl, n:int, newValue:varDbl.VarDbl, newVariance:varDbl.VarDbl,
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

    def __str__(self) -> str:
        return f'NotReliableException: {self.name} for {self.input} at {self.n}'

class NotStableException (Exception):
    def __init__(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[varDbl.VarDbl], inPrec:bool, outPrec:bool,
                 value:varDbl.VarDbl, variance:varDbl.VarDbl, n:int, newValue:varDbl.VarDbl, newVariance:varDbl.VarDbl,
                 *args: object) -> None:
        super().__init__(*args)
        self.input = input
        self.name = name
        self.s1dTaylor = s1dTaylor
        self.inPrec = inPrec

        self.value = value
        self.variance = variance
        self.n = n
        self.newValue = newValue
        self.newVariance = newVariance

    def __str__(self) -> str:
        return f'NotStableException: {self.name} for {self.input} at {self.n}'

class NotMonotonicException (Exception):
    def __init__(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[varDbl.VarDbl], inPrec:bool, outPrec:bool,
                 value:varDbl.VarDbl, variance:varDbl.VarDbl, n:int, newValue:varDbl.VarDbl, newVariance:varDbl.VarDbl,
                 *args: object) -> None:
        super().__init__(*args)
        self.input = input
        self.name = name
        self.s1dTaylor = s1dTaylor
        self.inPrec = inPrec

        self.value = value
        self.variance = variance
        self.n = n
        self.newValue = newValue
        self.newVariance = newVariance

    def __str__(self) -> str:
        return f'NotMonotonicException: {self.name} for {self.input} at {self.n}'


class DivergentException (Exception):
    def __init__(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[varDbl.VarDbl], inPrec:bool, outPrec:bool,
                 value:varDbl.VarDbl, variance:varDbl.VarDbl, n:int, newValue:varDbl.VarDbl, newVariance:varDbl.VarDbl,
                 *args: object) -> None:
        super().__init__(*args)
        self.input = input
        self.name = name
        self.s1dTaylor = s1dTaylor
        self.inPrec = inPrec

        self.value = value
        self.variance = variance
        self.n = n
        self.newValue = newValue
        self.newVariance = newVariance

    def __str__(self) -> str:
        return f'DivergentException: {self.name} for {self.input} at {self.n}'


class Taylor:
    TAU = 7.18e-7   # The stability test 
    MIN_TERMINATE_ORDER = 20

    __slots__ = ('_variance_threshold', '_momentum')

    def __init__(self, uncertainty_precision_threshold:float=None) -> None:
        self._momentum = momentum.Momentum()
        if uncertainty_precision_threshold is None:
            self._variance_threshold = 1.0 / self._momentum._binding / self._momentum._binding
        else:
            self._variance_threshold = uncertainty_precision_threshold * uncertainty_precision_threshold

    def taylor1d(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[typing.Union[float, varDbl.VarDbl]], 
                 inPrec:bool, outPrec:bool,
                 dumpPath=None, enableStabilityTruncation=True):
        '''
        1d Taylor expansion.
        
        @see     The paper for Variance Arithmetic on Taylor expansion convergence.
        
        @param name          The name of the Taylor expansion, for exception logging.
        @param s1dTaylor     The Taylor expansion coefficent, with f(x) as s1dTaylor[0]. It should already contains /n!.
        @param inPrec        If to expand by input precision
        @param outPrec       If the variance result needs to be multiplied by s1dTaylor[0].
        @param dumpPath                      If to dump the expansion to a file
        @param enableStabilityTruncation     If to truncate when the expansion becomes stable.
        
        @return  The result of taylor expansion with this as input.
        
        @exception InitException         If any item in s1dTaylor is not finite.
        @exception DivergentException    If the result is not finite.
        @exception NotReliableException  If the uncertainty of the variance is too large for its value. 
        @exceptopm NotMonotonicException If the result variance does not decrease monotonically. 
        @exceptopm NotStableException    If after maximal order expansion, the expansion is still not stable.       
        '''
        fw = None
        if dumpPath:
            fw = open(dumpPath, "w")
            fw.write("name\tvalue\tuncertainty\tvariance\tinPrec\toutPrec\tenableStabilityTruncation\tBinding\tMaxOrder\n")
            fw.write(f"{name}\t{input.value()}\t{input.uncertainty()}\t{input.variance()}"
                     f"\t{inPrec}\t{outPrec}\t{enableStabilityTruncation}"
                     f"\t{self._momentum._binding}\t{self._momentum._maxOrder}\n")
            for n in range(len(s1dTaylor)):
                fw.write(f"{n}\t")
            fw.write("\n")
            for n in range(len(s1dTaylor)):
                fw.write(f"{s1dTaylor[n].value()}\t")
            fw.write("\n")
            for n in range(len(s1dTaylor)):
                fw.write(f"{s1dTaylor[n].uncertainty()}\t")
            fw.write("\n")
            fw.write("2n\tExponent Value\tExponent Variance\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty")
            fw.write("\tlimit\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\n")

        value = varDbl.VarDbl(1, 0) if outPrec else s1dTaylor[0]
        variance = varDbl.VarDbl()
        var = varDbl.VarDbl(input.variance())
        if inPrec:
            var *= varDbl.VarDbl( 1 /input.value() /input.value() )
        varn = varDbl.VarDbl(var)
        prevVariance = None
        for n in range(2, min(len(s1dTaylor), self._momentum._maxOrder*2), 2):
            newValue = varn * s1dTaylor[n] * self._momentum.factor(n)
            newVariance = varDbl.VarDbl()
            for j in range(1, n):
                newVariance += varn * s1dTaylor[j] * s1dTaylor[n - j] * self._momentum.factor(n)
            for j in range(2, n, 2):
                newVariance -= varn * s1dTaylor[j] * self._momentum.factor(j) * \
                               s1dTaylor[n - j] * self._momentum.factor(n - j)
            value += newValue
            variance += newVariance
            unc = variance.value() *(Taylor.TAU**2)
            if fw:
                fw.write(f"{n}\t{varn.value()}\t{varn.variance()}"
                         f"\t{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}"
                         f"\t{unc}\t{newValue.value()}\t{newValue.variance()}\t{newVariance.value()}\t{newVariance.uncertainty()}\n")

            if (not math.isfinite(value.value())) or (not math.isfinite(value.variance())) \
                    or (not math.isfinite(variance.value())) or (not math.isfinite(variance.variance())):
                if fw:
                    fw.write("DivergentException\n")
                    fw.close()
                raise DivergentException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance)
            if variance.variance() > variance.value() * self._variance_threshold:
                if fw:
                    fw.write("NotReliableException\n")
                    fw.close()
                raise NotReliableException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance)

            if Taylor.MIN_TERMINATE_ORDER and (n >= Taylor.MIN_TERMINATE_ORDER) \
                    and (abs(prevVariance.value()) + unc < abs(newVariance.value())):
                if fw:
                    fw.write("NotMonotonicException\n")
                    fw.close()
                raise NotMonotonicException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance)
            varn *= var
            if varn.value() == 0:
                break
            if enableStabilityTruncation and (Taylor.MIN_TERMINATE_ORDER <= n) \
                    and ((abs(newVariance.value()) < unc) or (abs(newValue.value()) < math.ulp(value.value()))):
                break
            prevVariance = newVariance

        if enableStabilityTruncation and (n >= self._momentum._maxOrder*2):
            if fw:
                fw.write("NotStableException\n")
                fw.close()
            raise NotStableException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance)
        if outPrec:
            value *= s1dTaylor[0]
            variance *= s1dTaylor[0] * s1dTaylor[0]
        if fw:
            fw.write("Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\n")
            fw.write(f"{value.value()}\t{value.variance()}\t{variance.value()}\t{variance.variance()}\n")
            fw.close()
        return varDbl.VarDbl(value.value(), variance.value() + value.variance(), True)
    

    def polynominal(self, input:varDbl.VarDbl, sCoeff:tuple[typing.Union[float, varDbl.VarDbl]],
                    dumpPath=None):
        '''
        1d Taylor expansion for polynominal at "input" with "sCoeff".
        Allow input.value() +- input.uncertainty() to include 0
        '''
        if len(sCoeff) >= self._momentum._maxOrder:
            raise ValueError(f'The lenght {len(sCoeff)} of polynominal coefficient is more than {self._momentum._maxOrder}: {sCoeff}')
        exp = len(sCoeff) - 1
        s1dTaylor = [varDbl.VarDbl() if i else varDbl.VarDbl(sCoeff[0]) for i in range(2*exp + 1)]
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
        fw = None
        if dumpPath:
            fw = open(dumpPath, "w")
            fw.write("polynominal\tvalue\tuncertainty\tvariance\tBinding\tMaxOrder\n")
            fw.write(f"{len(sCoeff)}\t{input.value()}\t{input.uncertainty()}\t{input.variance()}"
                     f"\t{self._momentum._binding}\t{self._momentum._maxOrder}\nIndex\t")
            for n in range(len(s1dTaylor)):
                fw.write(f"{n}\t")
            fw.write("\nCoeff Value\t")
            for n in range(len(sCoeff)):
                if type(sCoeff[n]) == varDbl.VarDbl:
                    fw.write(f"{sCoeff[n].value()}\t")
                else:
                    fw.write(f"{sCoeff[n]}\t")
            fw.write("\nCoeff Uncertainty\t")
            for n in range(len(sCoeff)):
                if type(sCoeff[n]) == varDbl.VarDbl:
                    fw.write(f"{sCoeff[n].uncertainty()}\t")
                else:
                    fw.write("0\t")
            fw.write("\nTaylor Value\t")
            for n in range(len(s1dTaylor)):
                fw.write(f"{s1dTaylor[n].value()}\t")
            fw.write("\nTaylor Uncertainty\t")
            for n in range(len(s1dTaylor)):
                fw.write(f"{s1dTaylor[n].uncertainty()}\t")
            fw.write("\n")
            fw.write("2n\tExponent Value\tExponent Variance\tValue Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty")
            fw.write("\tNew Value Value\tNew Value Uncertainty\tNew Variance Value\tNew Variance Uncertainty\n")

        value = s1dTaylor[0] if type(s1dTaylor[0]) == varDbl.VarDbl else varDbl.VarDbl(s1dTaylor[0])
        variance = varDbl.VarDbl()
        var = varDbl.VarDbl(input.variance())
        varn = varDbl.VarDbl(var)
        for n in range(2, 2*exp + 1, 2):
            newValue = varn * s1dTaylor[n] * self._momentum.factor(n)
            newVariance = varDbl.VarDbl()
            for j in range(1, n):
                newVariance += varn * s1dTaylor[j] * s1dTaylor[n - j] * self._momentum.factor(n)
            for j in range(2, n, 2):
                newVariance -= varn * s1dTaylor[j] * self._momentum.factor(j) * \
                               s1dTaylor[n - j] * self._momentum.factor(n - j)
            value += newValue
            variance += newVariance
            if fw:
                fw.write(f"{n}\t{varn.value()}\t{varn.variance()}"
                         f"\t{value.value()}\t{value.uncertainty()}\t{variance.value()}\t{variance.uncertainty()}"
                         f"\t{newValue.value()}\t{newValue.variance()}\t{newVariance.value()}\t{newVariance.uncertainty()}\n")
            if (not math.isfinite(variance.value())) or (not math.isfinite(variance.variance())):
                if fw:
                    fw.write("DivergentException\n")
                    fw.close()
                raise varDbl.InitException(value, variance)
            varn *= var
            if varn.value() == 0:
                break
        if variance.variance() > variance.value() * self._variance_threshold:
            if fw:
                fw.write("NotStableException\n")
                fw.close()
            raise NotReliableException(input, f'{input}^{exp}', s1dTaylor, False, False,
                    value, variance, n, newValue, newVariance)
        if fw:
            fw.write("Value Value\tValue Uncertainty\tVariance Value\tVariance Uncertainty\n")
            fw.write(f"{value.value()}\t{value.variance()}\t{variance.value()}\t{variance.variance()}\n")
            fw.close()
        return varDbl.VarDbl(value.value(), variance.value() + value.variance(), True) 
    

    def exp(self) -> list[varDbl.VarDbl]:
        '''
        Generate s1dTaylor for taylor1d without s1dTaylor[0] for math.exp(var)
        '''
        sTaylor = []
        sTaylor.append(None)
        n = varDbl.VarDbl(1, 0)
        for i in range(1, self._momentum._maxOrder):
            n *= 1.0/i
            sTaylor.append(n)
        return sTaylor
     
    def log(self) -> list[varDbl.VarDbl]:
        sTaylor = []
        sTaylor.append(None)
        for i in range(1, self._momentum._maxOrder):
            sTaylor.append( varDbl.VarDbl((1/i) if ((i%2) == 1) else -1/i))
        return sTaylor


    def sin(self, x:float) -> list[float]:
        sTaylor = []
        sTaylor.append( varDbl.VarDbl(math.sin(x)) )
        n = varDbl.VarDbl(1, 0)
        for i in range(1, self._momentum._maxOrder):
            n *= 1/i
            match i % 4:
                case 0:
                    sTaylor.append(varDbl.VarDbl(math.sin(x) * n))
                case 1:
                    sTaylor.append(varDbl.VarDbl(math.cos(x) * n))  
                case 2:
                    sTaylor.append(varDbl.VarDbl(-math.sin(x) * n))
                case 3:
                    sTaylor.append(varDbl.VarDbl(-math.cos(x) * n))
        return sTaylor
    

    def power(self, exponent:float) -> list[float]:
        sTaylor = [varDbl.VarDbl(0, 0), varDbl.VarDbl(exponent)]
        for i in range(2, self._momentum._maxOrder):
            if sTaylor[-1].variance():
                sTaylor.append( sTaylor[-1] * ((exponent + 1)/i - 1) )
            else:
                sTaylor.append( varDbl.VarDbl(sTaylor[-1].value()/i * (exponent + 1 - i)) )
            if not sTaylor[-1]:
                break
        return sTaylor
    

