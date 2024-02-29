
import math

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

        self.value = value
        self.variance = variance
        self.n = n
        self.newValue = newValue
        self.newVariance = newVariance


class Taylor:
    __slots__ = ('_variance_threshold', '_momentum')

    def __init__(self, uncertainty_precision_threshold:float=None) -> None:
        self._momentum = momentum.Momentum()
        if uncertainty_precision_threshold is None:
            self._variance_threshold = 1.0 / self._momentum._binding / self._momentum._binding
        else:
            self._variance_threshold = uncertainty_precision_threshold * uncertainty_precision_threshold

    def taylor1d(self, input:varDbl.VarDbl, name:str, s1dTaylor:list[varDbl.VarDbl], inPrec:bool, outPrec:bool):
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
        for n in range(2, len(s1dTaylor), 2):
            newValue = s1dTaylor[n] * self._momentum.factor(n)
            newVariance = varDbl.VarDbl()
            for j in range(1, n):
                newVariance += s1dTaylor[j] * s1dTaylor[n - j] * self._momentum.factor(n)
            for j in range(2, n, 2):
                newVariance -= s1dTaylor[j] * self._momentum.factor(j) * \
                               s1dTaylor[n - j] * self._momentum.factor(n - j)
            value += newValue * varn
            variance += newVariance * varn
            if variance.variance() > variance.value() * self._variance_threshold:
                raise LossUncertaintyException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance)
            varn *= var
            if varn.value() == 0:
                break
        if outPrec:
            value *= s1dTaylor[0]
            variance *= s1dTaylor[0] * s1dTaylor[0]
        return varDbl.VarDbl(value.value(), variance.value() + value.variance(), True)

    def power1d(self, input:varDbl.VarDbl, exp:int, s1dTaylor:list[varDbl.VarDbl]):
        '''
        1d Taylor expansion for positive integer power.
        @return:            The output varDbl.VarDbl using Taylor expansion
        
        @param input:       The input varDbl.VarDbl    
        @param s1dTaylor:   The Taylor expansion coefficent, with f(x) as s1dTaylor[0].  
                            It should already contains /n!.

        Without any otptimization, or judgement of the convergence of the result.           
        '''
        def pow(p:int):
            return 1 if p <= 0 else math.pow(input.value(), p)
        value = s1dTaylor[0]
        variance = varDbl.VarDbl()
        var = varDbl.VarDbl(input.variance())
        varn = varDbl.VarDbl(var)
        for n in range(2, 2*exp + 1, 2):
            newValue = s1dTaylor[n] * self._momentum.factor(n) * pow(exp - n)
            newVariance = varDbl.VarDbl()
            for j in range(1, n):
                newVariance += s1dTaylor[j] * s1dTaylor[n - j] * self._momentum.factor(n)
            newVariance *= pow(2*exp - n)
            for j in range(2, n, 2):
                newVariance -= s1dTaylor[j] * self._momentum.factor(j) * pow(exp - j) * \
                               s1dTaylor[n - j] * self._momentum.factor(n - j) *pow(exp - (n - j))
            value += newValue * varn
            variance += newVariance * varn
            if variance.variance() > variance.value() * self._variance_threshold:
                raise LossUncertaintyException(input, f'{input}^{exp}', s1dTaylor, False, False,
                        value, variance, n, newValue, newVariance)
            varn *= var
            if varn.value() == 0:
                break
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
        sTaylor = []
        sTaylor.append(varDbl.VarDbl(0, 0))
        sTaylor.append(varDbl.VarDbl(exponent))
        for i in range(2, self._momentum._maxOrder):
            sTaylor.append( sTaylor[-1] * ((exponent + 1)/i - 1) )
        return sTaylor
    

