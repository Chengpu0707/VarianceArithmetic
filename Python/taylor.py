
import math
import typing

import varDbl
import momentum

class LossUncertaintyException (Exception):
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


class Taylor:
    TAU = 7.18e-7   # The stability test 

    __slots__ = ('_variance_threshold', '_momentum')

    def __init__(self, uncertainty_precision_threshold:float=None) -> None:
        self._momentum = momentum.Momentum()
        if uncertainty_precision_threshold is None:
            self._variance_threshold = 1.0 / self._momentum._binding / self._momentum._binding
        else:
            self._variance_threshold = uncertainty_precision_threshold * uncertainty_precision_threshold

    def taylor1d(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[varDbl.VarDbl], inPrec:bool, outPrec:bool,
                 enableStabilityTruncation=True):
        '''
        1d Taylor expansion.
        @return:            The output varDbl.VarDbl using Taylor expansion
        
        @param input:       The input varDbl.VarDbl    
        @param name:        The name of the Taylor expansion, for exception logging.
        @param s1dTaylor:   The Taylor expansion coefficent, with f(x) as s1dTaylor[0].  
                            It should already contains /n!.
        @param inPrec:      If to expand by input precision
        @param outPrec:     if the variance result needs to be multiplied by s1dTaylor[0]

        Without any otptimization, or judgement of the convergence of the result.           
        '''
        value = varDbl.VarDbl(1, 0) if outPrec else s1dTaylor[0]
        variance = varDbl.VarDbl()
        var = varDbl.VarDbl(input.variance())
        if inPrec:
            var *= varDbl.VarDbl( 1 /input.value() /input.value() )
        varn = varDbl.VarDbl(var)
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
            if variance.variance() > variance.value() * self._variance_threshold:
                raise LossUncertaintyException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance)
            varn *= var
            if varn.value() == 0:
                break
            if enableStabilityTruncation:
                unc = variance.value()*Taylor.TAU
                if (math.sqrt(abs(newVariance.value())) < unc) and \
                        (abs(newValue.value()) < max(unc, math.ulp(value.value()))):
                    break
        if enableStabilityTruncation and (n >= self._momentum._maxOrder*2):
            raise NotStableException(input, name, s1dTaylor, inPrec, outPrec,
                    value, variance, n, newValue, newVariance)
        if outPrec:
            value *= s1dTaylor[0]
            variance *= s1dTaylor[0] * s1dTaylor[0]
        return varDbl.VarDbl(value.value(), variance.value() + value.variance(), True)

    def polynominal(self, input:varDbl.VarDbl, sCoeff:tuple[typing.Union[float, varDbl.VarDbl]]):
        '''
        1d Taylor expansion for polynominal at "input" with "sCoeff".
        Allow input.value() +- input.uncertainty() to include 0
        '''
        if len(sCoeff) >= self._momentum._maxOrder:
            raise ValueError(f'The lenght {len(sCoeff)} of polynominal coefficient is more than {self._momentum._maxOrder}: {sCoeff}')
        exp = len(sCoeff) - 1
        s1dTaylor = [varDbl.VarDbl() for i in range(2*exp + 1)]
        sPow:map[int, float] = {1: input.value()}

        def pow(p:int):
            p = int(p)
            if p <= 0:
                return 1
            if (p not in sPow):
                maxP = max(sPow)
                while (maxP < p):
                    sPow[maxP + 1] = sPow[maxP] * input.value()
                    maxP += 1
            return sPow[p]

        for j, coeff in enumerate(sCoeff):
            if not coeff:
                continue
            sTaylor = [1, j]
            for k in range(2, j + 1):
                sTaylor.append( sTaylor[-1] * (j + 1 - k)//k )
            for k in range(j + 1):
                s1dTaylor[k] += coeff * sTaylor[k] * pow(j - k)
        value = s1dTaylor[0]
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
            if (not math.isfinite(variance.value())) or (not math.isfinite(variance.variance())):
                raise varDbl.UncertaintyException(value, variance)
            varn *= var
            if varn.value() == 0:
                break
        if variance.variance() > variance.value() * self._variance_threshold:
            raise LossUncertaintyException(input, f'{input}^{exp}', s1dTaylor, False, False,
                    value, variance, n, newValue, newVariance)
        return varDbl.VarDbl(value.value(), variance.value() + value.variance(), True) if type(value) == varDbl.VarDbl \
                    else varDbl.VarDbl(value, variance.value(), True)
    

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
            sTaylor.append( varDbl.VarDbl(1 if ((i%2) == 1) else -1, 0) * (1/i) )
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
            sTaylor.append( sTaylor[-1] * ((exponent + 1)/i - 1) )
            if not sTaylor[-1]:
                break
        return sTaylor
    

