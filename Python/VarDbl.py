import math
from typing import Optional, Union

class ValDblValueError (BaseException):
    def __init__(self, value: float, *args: object) -> None:
        super().__init__(*args)
        self.value = value

class ValDblUncertaintyError (BaseException):
    def __init__(self, value: float, uncertainty: float, *args: object) -> None:
        super().__init__(*args)
        self.value = value
        self.uncertainty = uncertainty


class VarDbl:
    __slots__ = ('_value', '_variance')

    def value(self):
        return self._value
    
    def uncertainty(self):
        return math.sqrt(self._variance)
        
    def __init__(self, value: Union[float, int], uncertainty: Optional[float]=None) -> None:
        '''
        Intialize with "value" and "uncertainty".
        "uncertainty" will be absolute, and limited between 
            math.sqrt(sys.float_info.min) and math.sqrt(sys.float_info.max)
        Both have to be finite. 
            Otherwise ValDblValueError or ValDblUncertaintyError will throw.

        When "uncertainty" is not specified:
             *) An int is initialized with uncertainty = 0
             *) An float is initialized with uncertainty = math.ulp
        '''
        if uncertainty is None:
            if type(value) == int:
                if value == int(float(value)):
                    self._value = float(value)
                    self._variance = 0.0
                    return
                value = float(value)
            uncertainty = math.ulp(value) 
        if not math.isfinite(value):
            raise ValDblValueError(value)
        variance = uncertainty * uncertainty
        if not math.isfinite(variance):
            raise ValDblUncertaintyError(value, uncertainty)
        self._value = value
        self._variance = variance

    def __str__(self) -> str:
        return f'{self.value():.6e}~{self.uncertainty():.3e}'

    def __repr__(self) -> str:
        return f'{repr(self._value)}~{repr(self._variance)}'

