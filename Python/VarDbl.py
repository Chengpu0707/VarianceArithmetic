import math
import typing


class ValueException (Exception):
    def __init__(self, value: float, *args: object) -> None:
        super().__init__(*args)
        self.value = value

class UncertaintyException (Exception):
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
    DEVIATION_OF_LSB = 1.0 / math.sqrt(3)
        # rounding error is uniformly distrubuted within LSB of float
    # Taylor expansion parameters are defined in taylor.Taylor and momentum.Momentum
    _taylor = None

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
        If "bUncertaintyAsVariance" is Ture, the uncertainty actually means variance, 
            which should only be True during intermediate calculations. 
        Both value and variance have to be finite. 
            Otherwise ValueException or UncertaintyException will throw.

        When "uncertainty" is not specified:
             *) An int is initialized with uncertainty = 0
             *) An float is initialized with uncertainty = math.ulp
        '''
        if uncertainty is None:
            if type(value) == VarDbl:
                self._value = value._value
                self._variance = value._variance
                return
            if type(value) == int:
                if value == int(float(value)):
                    self._value = float(value)
                    self._variance = 0.0
                    return
                value = float(value)
            uncertainty = math.ulp(value) * VarDbl.DEVIATION_OF_LSB
        variance = uncertainty if bUncertaintyAsVariance else uncertainty * uncertainty
        if not math.isfinite(value):
            raise ValueException(value, "__init__")
        if not math.isfinite(variance):
            raise UncertaintyException(value, uncertainty, "__init__")
        self._value = float(value)
        self._variance = float(variance)

    def __str__(self) -> str:
        return f'{self.value():.6e}~{self.uncertainty():.3e}'

    def __repr__(self) -> str:
        return f'{repr(self.value())}~{repr(self.variance())}'
    

    def __add__(self, other):
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        value = self.value() + other.value()
        variance = self.variance() + other.variance()
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
        if VarDbl._taylor is None:
            import taylor
            VarDbl._taylor = taylor.Taylor()  
        s1dTaylor = VarDbl._taylor.exp()
        s1dTaylor[0] = math.exp(var.value())
        return VarDbl._taylor.taylor1d(var, "exp", s1dTaylor, False, True)
    
    @staticmethod
    def log(var):
        if VarDbl._taylor is None:
            import taylor
            VarDbl._taylor = taylor.Taylor()  
        s1dTaylor = VarDbl._taylor.log()
        s1dTaylor[0] = math.log(var.value())
        return VarDbl._taylor.taylor1d(var, f"log({var})", s1dTaylor, True, False)
    
    @staticmethod
    def sin(var):
        if VarDbl._taylor is None:
            import taylor
            VarDbl._taylor = taylor.Taylor()  
        s1dTaylor = VarDbl._taylor.sin(var.value())
        s1dTaylor[0] = math.sin(var.value())
        return VarDbl._taylor.taylor1d(var, f"sin({var})", s1dTaylor, False, False)
    
    def __pow__(self, exp):
        match exp:
            case 0:
                return VarDbl(1, 0)
            case 1:
                return VarDbl(self)

        if VarDbl._taylor is None:
            import taylor
            VarDbl._taylor = taylor.Taylor()  
        exp = float(exp)
        s1dTaylor = VarDbl._taylor.power(exp)
        s1dTaylor[0] = math.pow(self.value(), exp)
        if (exp > 0) and (math.ceil(exp) == math.floor(exp)):
            return VarDbl._taylor.power1d(self, int(exp), s1dTaylor)
        return VarDbl._taylor.taylor1d(self, f"{self}**{exp}", s1dTaylor, True, True)
    
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
        deltaValue = delta=math.ulp(var.value())
    self.assertAlmostEqual(value, var.value(), delta=deltaValue)
    if uncertainty is None:
         uncertainty = math.ulp(var.value())
    if deltaUncertainty is None:
         deltaUncertainty = math.ulp(max(uncertainty, var.uncertainty()))
    self.assertAlmostEqual(uncertainty, var.uncertainty(), delta=deltaUncertainty)
