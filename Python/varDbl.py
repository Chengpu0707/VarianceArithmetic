import math
import typing

'''
VarDbl uses Taylor which import VarDbl, causing a circular dependency, 
because of ** operator, and VarDbl.<func>().

The current solution is to import Taylor within VarDbl.taylor()

'''

class InitException (Exception):
    def __init__(self, value: float, uncertainty: float, *args: object) -> None:
        super().__init__(*args)
        self.value = value
        self.uncertainty = uncertainty


class VarDbl:
    '''
    A class for variance arithmetic that contains basic operations.

    The Taylor expansion part is in the class Taylor
    '''
    BINDING_FOR_EQUAL = 0.67448975
        # z value for 50% probability of equal
    
    DOUBLE_MAX_SIGNIFICAND = (1 << 53) - 1
    DOUBLE_MAX_PRECISE_SIGNIFICAND_BIT = 30
    DEVIATION_OF_LSB = 1.0 / math.sqrt(3)
        # rounding error is uniformly distrubuted within LSB of float
    @staticmethod
    def ulp(value:typing.Union[float, int]) -> float:
        if type(value) == float:
            return math.ulp(value) * VarDbl.DEVIATION_OF_LSB
        if type(value) == int:
            round = 0.0
            posi = True
            val = abs(value)
            while VarDbl.DOUBLE_MAX_SIGNIFICAND < val:
                if val & 1:
                    if posi:
                        round += 1
                        posi = False
                    else:
                        round -= 1
                        posi = True
                round *= 0.5
                val >>= 1
            return round
        return VarDbl.ulp(float(value))
            

    # Taylor expansion parameters are defined in taylor.Taylor and momentum.Momentum
    _taylor = None
    @staticmethod
    def taylor():
        if VarDbl._taylor is None:
            import taylor
            VarDbl._taylor = taylor.Taylor()
        return VarDbl._taylor


    __slots__ = ('_value', '_variance')

    def value(self):
        return self._value
    
    def uncertainty(self):
        return math.sqrt(self.variance())
    
    def variance(self):
        return self._variance
        
    def __init__(self, value: typing.Union[float, int]=0, 
                 uncertainty: typing.Optional[float]=None,
                 bUncertaintyAsVariance=False) -> None:
        '''
        Intialize with "value" and "uncertainty".
        "uncertainty" will be absolute, and limited between 
            math.sqrt(sys.float_info.min) and math.sqrt(sys.float_info.max)
        If "bUncertaintyAsVariance" is True, the uncertainty actually means variance, 
            which should only be True during intermediate calculations. 
        Both value and variance have to be finite. 
            Otherwise InitException will throw.

        When "uncertainty" is not specified:
             *) An int which is not more than DOUBLE_MAX_SIGNIFICAND, is initialized with uncertainty = 0.
             *) An int which is more than DOUBLE_MAX_SIGNIFICAND, is initialized with uncertainty = rounding error.
             *) A float which is 2's fraction larger than 2^{-40}, is initialized with uncertainty = 0.
             *) Otherwise, a float is initialized with uncertainty = math.ulp.
        '''
        if uncertainty is None:
            if type(value) == VarDbl:
                self._value = value._value
                self._variance = value._variance
                return
            if type(value) == int:
                uncertainty = VarDbl.ulp(value)
            else:
                value = float(value)
                if not math.isfinite(value):
                    raise InitException(value, None)
                for bit in (0, VarDbl.DOUBLE_MAX_PRECISE_SIGNIFICAND_BIT + 1):
                    val = value * (1 << (VarDbl.DOUBLE_MAX_PRECISE_SIGNIFICAND_BIT - bit))
                    if math.isfinite(val):
                        break
                uncertainty = 0 if math.floor(val) == math.ceil(val) else VarDbl.ulp(value)
        variance = uncertainty if bUncertaintyAsVariance else uncertainty * uncertainty
        if (not math.isfinite(value)) or (not math.isfinite(variance)):
            raise InitException(value, uncertainty)
        self._value = float(value)
        self._variance = float(variance)

    def __str__(self) -> str:
        return f'{self.value():.6e}~{self.uncertainty():.3e}'

    def __repr__(self) -> str:
        return f'{repr(self.value())}~{repr(self.variance())}'
    
    def __bool__(self) -> bool:
        return (self._value != 0) or (self._variance != 0)

    def __add__(self, other):
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        value = self.value() + other.value()
        variance = self.variance() + other.variance()
        if (not math.isfinite(value)) or (not math.isfinite(variance)):
            raise InitException(f"{self} + {other} = {value}~{variance}", value, math.sqrt(variance))
        if (not variance) and (abs(self.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) \
                          and (abs(other.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) \
                          and (VarDbl.DOUBLE_MAX_SIGNIFICAND <= abs(value)):
            return VarDbl(int(self.value()) + int(other.value()))
        return VarDbl(value, variance, True)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        ret = VarDbl()
        ret._value = - self._value
        ret._variance = self._variance
        return ret

    def __sub__(self, other):
        return self + other.__neg__()
    
    def __rsub__(self, other):
        return -(self - other)
    

    def __mul__(self, other):
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        value = self.value() * other.value()
        variance = self.variance() * other.value() * other.value() +\
                    other.variance() *self.value() * self.value() +\
                    self.variance() * other.variance()
        if (not math.isfinite(value)) or (not math.isfinite(variance)):
            raise InitException(f"{self} * {other} = {value}~{variance}", value, math.sqrt(variance))
        if (not variance) and (abs(self.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) \
                          and (abs(other.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) \
                          and (VarDbl.DOUBLE_MAX_SIGNIFICAND <= abs(value)):
            return VarDbl(int(self.value()) * int(other.value()))
        return VarDbl(value, variance, True)
    
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other: object):
        if type(other) != VarDbl:
            other = VarDbl(other)
        return self * (other ** -1)

    def __rtruediv__(self, other: object):
        other = VarDbl(value=other)
        return other / self

    @staticmethod
    def exp(var):
        s1dTaylor = VarDbl.taylor().exp()
        s1dTaylor[0] = math.exp(var.value())
        return VarDbl.taylor().taylor1d(var, "exp", s1dTaylor, False, True)
    
    @staticmethod
    def log(var):
        s1dTaylor = VarDbl.taylor().log()
        s1dTaylor[0] = math.log(var.value())
        return VarDbl.taylor().taylor1d(var, f"log({var})", s1dTaylor, True, False)
    
    @staticmethod
    def sin(var):
        s1dTaylor = VarDbl.taylor().sin(var.value())
        s1dTaylor[0] = math.sin(var.value())
        return VarDbl.taylor().taylor1d(var, f"sin({var})", s1dTaylor, False, False)
    
    def __pow__(self, exp):
        match exp:
            case 0:
                return VarDbl(1, 0)
            case 1:
                return VarDbl(self)
        if (exp > 0) and ((type(exp) == int) or (math.ceil(exp) == math.floor(exp))):
            sCoeff = [0] * int(exp)
            sCoeff.append(1)
            return VarDbl.taylor().polynominal(self, sCoeff)
        s1dTaylor = VarDbl.taylor().power(exp)
        s1dTaylor[0] = math.pow(self.value(), exp)
        return VarDbl.taylor().taylor1d(self, f"{self}**{exp}", s1dTaylor, True, True)
    
    def __eq__(self, other: object) -> bool:
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        diff = self - other
        if diff.value() == 0:
            return True
        uncertainty = diff.uncertainty()
        if uncertainty == 0:
            return False
        return abs(diff.value()/uncertainty) <= VarDbl.BINDING_FOR_EQUAL

    def __lt__(self, other: object) -> bool:
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        if self == other:
            return False
        return self.value() < other.value()

    def __gt__(self, other: object) -> bool:
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        if self == other:
            return False
        return self.value() > other.value()

    def __le__(self, other: object) -> bool:
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        if self == other:
            return True
        return self.value() < other.value()

    def __ge__(self, other: object) -> bool:
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        if self == other:
            return True
        return self.value() > other.value()
  







def validate(self, var:VarDbl, value:float, uncertainty:float=None,
             deltaValue=None, deltaUncertainty=None):
    '''
    "self" should refer to a unittest.TestCase instance
    '''
    if deltaValue is None:
        deltaValue = math.ulp(var.value())
    self.assertAlmostEqual(value, var.value(), delta=deltaValue)
    if uncertainty is None:
         uncertainty = math.ulp(var.value())
    if deltaUncertainty is None:
         deltaUncertainty = math.ulp(max(uncertainty, var.uncertainty()))
    self.assertAlmostEqual(uncertainty, var.uncertainty(), delta=deltaUncertainty)
