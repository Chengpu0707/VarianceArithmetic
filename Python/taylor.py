
import math

from varDbl import VarDbl
from momentum import Momentum

class LossUncertaintyException (Exception):
    def __init__(self, input:VarDbl, name:str, s1dTaylor:list[VarDbl], inPrec:bool, outPrec:bool,
                 value:VarDbl, variance:VarDbl, n:int, newValue:VarDbl, newVariance:VarDbl,
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
        self._momentum = Momentum()
        if uncertainty_precision_threshold is None:
            self._variance_threshold = 1.0 / self._momentum._binding / self._momentum._binding
        else:
            self._variance_threshold = uncertainty_precision_threshold * uncertainty_precision_threshold

    def taylor1d(self, input:VarDbl, name:str, s1dTaylor:list[VarDbl], inPrec:bool, outPrec:bool):
        '''
        1d Taylor expansion.
        @return:            The output VarDbl using Taylor expansion
        
        @param input:       The input VarDbl    
        @param name:        The name of the Taylor expansion, for exception logging.
        @param s1dTaylor:   The Taylor expansion coefficent, with f(x) as s1dTaylor[0].  
                            It should already contains /n!.
        @param inPrec:      If to expand by input precision
        @param outPrec:     if the variance result needs to be multiplied by s1dTaylor[0]

        Without any otptimization, or judgement of the convergence of the result.           
        '''
        value = VarDbl(1, 0) if outPrec else s1dTaylor[0]
        variance = VarDbl()
        var = VarDbl(input.variance())
        if inPrec:
            var *= VarDbl( 1 /input.value() /input.value() )
        varn = VarDbl(var)
        for n in range(2, self._momentum._maxOrder, 2):
            newValue = s1dTaylor[n] * self._momentum.factor(n) * varn
            newVariance = VarDbl()
            for j in range(1, n):
                newVariance += s1dTaylor[j] * s1dTaylor[n - j] * self._momentum.factor(n) * varn
            for j in range(2, n, 2):
                newVariance -= s1dTaylor[j] * self._momentum.factor(j) * \
                            s1dTaylor[n - j] * self._momentum.factor(n - j) * \
                            varn
            value += newValue
            variance += newVariance
            if variance.variance() > variance.value() * self._variance_threshold:
                raise LossUncertaintyException(input, name, s1dTaylor, inPrec, outPrec,
                        value, variance, n, newValue, newVariance)
            varn *= var
        if outPrec:
            value *= s1dTaylor[0]
            variance *= s1dTaylor[0] * s1dTaylor[0]
        return VarDbl(value.value(), variance.value() + value.variance(), True)

    
    def exp(self) -> list[VarDbl]:
        '''
        Generate s1dTaylor for taylor1d without s1dTaylor[0] for math.exp(var)
        '''
        sTaylor = []
        sTaylor.append(None)
        n = VarDbl(1, 0)
        for i in range(1, self._momentum._maxOrder):
            n *= 1.0/i
            sTaylor.append(n)
        return sTaylor
     
    def log(self) -> list[VarDbl]:
        sTaylor = []
        sTaylor.append(None)
        for i in range(1, self._momentum._maxOrder):
            sTaylor.append( VarDbl(1 if ((i%2) == 1) else -1, 0) * (1/i) )
        return sTaylor


    def sin(self, x:float) -> list[float]:
        sTaylor = []
        sTaylor.append( VarDbl(math.sin(x)) )
        n = VarDbl(1, 0)
        for i in range(1, self._momentum._maxOrder):
            n *= 1/i
            match i % 4:
                case 0:
                    sTaylor.append(VarDbl(math.sin(x) * n))
                case 1:
                    sTaylor.append(VarDbl(math.cos(x) * n))  
                case 2:
                    sTaylor.append(VarDbl(-math.sin(x) * n))
                case 3:
                    sTaylor.append(VarDbl(-math.cos(x) * n))
        return sTaylor
    
    def power(self, exponent:float) -> list[float]:
        sTaylor = []
        sTaylor.append(VarDbl(0, 0))
        sTaylor.append(VarDbl(exponent))
        exponent -= 1
        for i in range(2, self._momentum._maxOrder):
            sTaylor.append( sTaylor[-1] * (exponent / i) )
            exponent -= 1
        return sTaylor
    

