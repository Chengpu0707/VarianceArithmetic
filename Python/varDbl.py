import math
import numbers
import sys
import typing
import unittest

'''
VarDbl uses Taylor which import VarDbl, causing a circular dependency, 
because of ** operator, and VarDbl.<func>().

The current solution is to import Taylor within VarDbl.taylor()

'''

class InitException (Exception):
    def __init__(self, value: float, variance: float, *args: object) -> None:
        super().__init__(*args)
        self.value = value
        self.variance = variance


class VarDbl (numbers.Number):
    '''
    A class for variance arithmetic which contains basic operations including Taylor expansion.
        The Taylor expansion part is in the class Taylor
    Each value represents a specified random variable containing:
        A value(): the expected mean 
        A uncertainty(): the expected deviation, with its square as variance()
    
    Two VarDbl objects are equal if the absolute z-value of their difference is less than {BINDING_FOR_EQUAL}.
        For this reason, self.__hash__() raise NotImplenentedException, so that a VarDbl can not be a key.
    '''
    BINDING_FOR_EQUAL = 0.67448975
        # z value for 50% probability of equal
    DOUBLE_MAX_PRECISE_FILTER = (1 << 23) - 1
    
    DOUBLE_MAX_SIGNIFICAND = (1 << 53) - 1
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
            
    __slots__ = ('_value', '_uncertainty')

    def value(self):
        return self._value
    
    def uncertainty(self):
        return self._uncertainty
    
    def variance(self):
        return self._uncertainty**2
    
    def precision(self):
        try:
            return self.uncertainty() / abs(self.value())
        except ZeroDivisionError:
            return float('inf')

        
    def __init__(self, value: typing.Union[float, int, str]=0, 
                 uncertainty: typing.Union[float, str, None]=None) -> None:
        '''
        Intialize with "value" and "uncertainty".
        "uncertainty" will be absolute, and limited between 
            math.sqrt(sys.float_info.min) and math.sqrt(sys.float_info.max)
        Both value and vuncertainty have to be finite. 
            Otherwise InitException will throw.

        When "uncertainty" is not specified:
             *) An int which is not more than DOUBLE_MAX_SIGNIFICAND, is initialized with uncertainty = 0.
             *) An int which is more than DOUBLE_MAX_SIGNIFICAND, is initialized with uncertainty = rounding error.
             *) A float which is 2's fraction larger than 2^{-40}, is initialized with uncertainty = 0.
             *) Otherwise, a float is initialized with uncertainty = math.ulp.
        '''
        if isinstance(value, str):
            value = float(value)
        if isinstance(uncertainty, str):
            uncertainty = float(uncertainty)
        if uncertainty is None:
            if type(value) == VarDbl:
                self._value = value._value
                self._uncertainty = value._uncertainty
                return
            if type(value) == int:
                uncertainty = VarDbl.ulp(value)
            else:
                value = float(value)
                if not math.isfinite(value):
                    raise InitException(value, None, 'Init value=inf')
                sig, _ = math.frexp(value)
                sig = int(sig * (VarDbl.DOUBLE_MAX_SIGNIFICAND + 1))
                uncertainty = VarDbl.ulp(value) if (sig & VarDbl.DOUBLE_MAX_PRECISE_FILTER) != 0 else 0
        else:
            uncertainty = abs(uncertainty)
        if (not math.isfinite(value)) or (not math.isfinite(uncertainty)):
            raise InitException(value, uncertainty, f'Init value={value}, uncertainty={uncertainty}')
        self._value = float(value)
        self._uncertainty = float(uncertainty)

    def __str__(self) -> str:
        return f'{self.value():.6e}~{self.uncertainty():.3e}'

    def __repr__(self) -> str:
        return f'{repr(self.value())}~{repr(self.uncertainty())}'
    
    def __bool__(self) -> bool:
        return (self._value != 0) or ( (self._uncertainty != 0))
    
    def __abs__(self):
        return VarDbl(abs(self.value()), self.uncertainty())
    
    def __add__(self, other):
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        value = self.value() + other.value()
        if self.uncertainty() == 0:
            uncertainty = other.uncertainty()
        elif other.uncertainty() == 0:
            uncertainty = self.uncertainty()
        else:
            try:
                uncertainty = math.sqrt(self.variance() + other.variance())
            except OverflowError:
                raise InitException(value, uncertainty, f"{self} + {other}")
        return VarDbl(value, uncertainty)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        ret = VarDbl()
        ret._value = - self._value
        ret._uncertainty = self._uncertainty
        return ret

    def __sub__(self, other):
        return self + other.__neg__()
    
    def __rsub__(self, other):
        return -(self - other)
    

    def __mul__(self, other):
        if type(other) != VarDbl:
            other = VarDbl(value=other)
        try:
            value = self.value() * other.value()
        except OverflowError:
            raise InitException(None, None, f"{self} * {other}")
        try:
            variance = self.variance() * other.value() * other.value() +\
                        other.variance() *self.value() * self.value() +\
                        self.variance() * other.variance()
        except OverflowError:
            raise InitException(value, None, f"{self} * {other}")
        if (not math.isfinite(value)) or (not math.isfinite(variance)):
            raise InitException(value, variance, f"{self} * {other} = {value}~{variance}")
        if (not variance) and (abs(self.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) \
                          and (abs(other.value()) < VarDbl.DOUBLE_MAX_SIGNIFICAND) \
                          and (VarDbl.DOUBLE_MAX_SIGNIFICAND <= abs(value)):
            return VarDbl(int(self.value()) * int(other.value()))
        return VarDbl(value, math.sqrt(variance))
    
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other: object):
        if type(other) != VarDbl:
            other = VarDbl(other)
        return self * (other ** -1)

    def __rtruediv__(self, other: object):
        other = VarDbl(value=other)
        return other / self
   
    def __pow__(self, exp):
        import taylor
        return taylor.Taylor.pow(self, exp)
    
    def __hash__(self) -> int:
        raise NotImplemented('Difficult to find hash')
    
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
  

SMALLEST_SUB_NORAMAL_FLOAT = sys.float_info.min * sys.float_info.epsilon

def assertVarDblEqual(self:unittest.TestCase, l:VarDbl, r:VarDbl, valPrec=1e-6, uncPrec=1e-6):
    try:
        if (r.value() == 0):
            self.assertLessEqual(abs(l.value()), SMALLEST_SUB_NORAMAL_FLOAT)
        elif (l.value() == 0):
            self.assertLessEqual(abs(r.value()), SMALLEST_SUB_NORAMAL_FLOAT)
        else:
            self.assertAlmostEqual(l.value() / r.value(), 1, delta=valPrec)
    except AssertionError as ex:
        raise ex
    try:
        if (r.uncertainty() == 0):
            self.assertLessEqual(abs(l.uncertainty()), SMALLEST_SUB_NORAMAL_FLOAT)
        elif (l.uncertainty() == 0):
            self.assertLessEqual(abs(r.uncertainty()), SMALLEST_SUB_NORAMAL_FLOAT)
        else:
            self.assertAlmostEqual(l.uncertainty() / r.uncertainty(), 1, delta=uncPrec)
    except AssertionError as ex:
        raise ex




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
