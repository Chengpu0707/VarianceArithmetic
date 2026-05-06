import enum
import functools
import itertools
import math
import operator
import sympy
import typing

import momentum


class EDistrType(enum.Enum):
    Uniform = enum.auto()
    Gaussian = enum.auto()


class zeta(sympy.Function):
    """Symbolic bounded moment ζ(n, κ) of the truncated Normal distribution."""
    nargs = 2

    @classmethod
    def eval(cls, n, _):
        if n.is_odd:
            return sympy.Integer(0)

    def _eval_evalf(self, _):
        n, kappa = self.args
        if n.is_integer and kappa.is_number:
            return sympy.Float(momentum.Normal(bounding=float(kappa))[int(n)])


class InVarException(Exception):
    pass


class InVar:
    __slots__ = ('_value', '_deviation', '_distr_type', '_kappa', '_samples')

    def __init__(self, value: sympy.Symbol, deviation: sympy.Symbol,
                 distr_type: EDistrType, kappa: float = None, samples: int = 10000):
        if not isinstance(value, sympy.Symbol):
            raise InVarException(f'value must be a sympy.Symbol, got {type(value)}')
        if not isinstance(deviation, sympy.Symbol):
            raise InVarException(f'deviation must be a sympy.Symbol, got {type(deviation)}')
        if not isinstance(distr_type, EDistrType):
            raise InVarException(f'distr_type must be an EDistrType, got {type(distr_type)}')
        if kappa is None:
            kappa = 5.0 if distr_type == EDistrType.Gaussian else math.sqrt(3)
        if not isinstance(kappa, float):
            raise InVarException(f'kappa must be a float, got {type(kappa)}')
        if kappa <= 0:
            raise InVarException(f'kappa must be positive, got {kappa}')
        if distr_type == EDistrType.Uniform and kappa > math.sqrt(3):
            raise InVarException(f'kappa must be <= sqrt(3) for Uniform distribution, got {kappa}')
        if not isinstance(samples, int) or samples <= 0:
            raise InVarException(f'samples must be a positive int, got {samples}')
        self._value = value
        self._deviation = deviation
        self._distr_type = distr_type
        self._kappa = kappa
        self._samples = samples

    @property
    def value(self) -> sympy.Symbol:
        return self._value

    @property
    def deviation(self) -> sympy.Symbol:
        return self._deviation

    @property
    def distr_type(self) -> EDistrType:
        return self._distr_type

    @property
    def kappa(self) -> float:
        return self._kappa

    @property
    def samples(self) -> int:
        return self._samples

    def moment(self, order: int) -> typing.Union[float, typing.Callable[[int, float], float]]:
        if self._distr_type == EDistrType.Gaussian:
            return zeta(order, self._kappa)
        if order % 2 == 1:
            return 0.0
        # Formula (2.22): ζ(n) = 2ρ(κ) · κ^(n+1) / (n+1), 2ρ(κ) = 1/√3
        return self._kappa ** (order + 1) / (math.sqrt(3) * (order + 1))


def _multi_indices(n: int, k: int):
    """Yield all length-n tuples of non-negative ints summing to k."""
    if n == 1:
        yield (k,)
        return
    for i in range(k + 1):
        for rest in _multi_indices(n - 1, k - i):
            yield (i,) + rest


def _multi_indices_with_min1(n_vars: int, total: int):
    """Yield length-n_vars tuples with each element >= 1 summing to total."""
    if total < n_vars:
        return
    for p in _multi_indices(n_vars, total - n_vars):
        yield tuple(pk + 1 for pk in p)


class TaylorException(Exception):
    pass


class Taylor:
    __slots__ = ('_function', '_in_vars', '_max_order', '_coeffs')

    def __init__(self, function: sympy.Expr, in_vars: tuple, max_order: int = 200):
        if not isinstance(function, sympy.Expr):
            raise TaylorException(f'function must be a sympy.Expr, got {type(function)}')
        if not isinstance(in_vars, tuple):
            raise TaylorException(f'in_vars must be a tuple, got {type(in_vars)}')
        if len(in_vars) == 0:
            raise TaylorException('in_vars must have at least one element')
        for i, v in enumerate(in_vars):
            if not isinstance(v, InVar):
                raise TaylorException(f'in_vars[{i}] must be an InVar, got {type(v)}')
        deviation_symbols = {v.deviation for v in in_vars}
        if function.free_symbols & deviation_symbols:
            raise TaylorException(
                f'function must not contain deviation symbols: {function.free_symbols & deviation_symbols}')
        if not isinstance(max_order, int) or isinstance(max_order, bool) or max_order < 0:
            raise TaylorException(f'max_order must be a non-negative int, got {max_order}')
        self._function = function
        self._in_vars = in_vars
        self._max_order = max_order
        symbols = [v.value for v in in_vars]
        n_vars = len(symbols)
        zero = (0,) * n_vars
        self._coeffs = {zero: function}
        for k in range(1, max_order + 1):
            for alpha in _multi_indices(n_vars, k):
                for i, ai in enumerate(alpha):
                    if ai > 0:
                        parent = alpha[:i] + (ai - 1,) + alpha[i + 1:]
                        self._coeffs[alpha] = sympy.diff(self._coeffs[parent], symbols[i]) / ai
                        break

    @property
    def function(self) -> sympy.Expr:
        return self._function

    @property
    def in_vars(self) -> tuple:
        return self._in_vars

    @property
    def max_order(self) -> int:
        return self._max_order

    @property
    def coeffs(self) -> dict:
        return self._coeffs

    def at(self, *orders: int) -> sympy.Expr:
        if len(orders) != len(self._in_vars):
            raise TaylorException(
                f'expected {len(self._in_vars)} orders, got {len(orders)}')
        for o in orders:
            if not isinstance(o, int) or isinstance(o, bool) or o < 0:
                raise TaylorException(f'each order must be a non-negative int, got {o}')
        total = sum(orders)
        if total > self._max_order:
            raise TaylorException(f'total order {total} exceeds max_order {self._max_order}')
        return self._coeffs[orders]

    def varAt(self, *orders: int) -> sympy.Expr:
        if len(orders) != len(self._in_vars):
            raise TaylorException(
                f'expected {len(self._in_vars)} orders, got {len(orders)}')
        for o in orders:
            if not isinstance(o, int) or isinstance(o, bool) or o < 1:
                raise TaylorException(f'each order must be a positive int, got {o}')
        total = sum(orders)
        if total > self._max_order:
            raise TaylorException(f'total order {total} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        p = orders
        dev_prod = functools.reduce(operator.mul,
                                    (self._in_vars[k].deviation ** p[k] for k in range(N)),
                                    sympy.Integer(1))
        full_moment = functools.reduce(operator.mul,
                                       (self._in_vars[k].moment(p[k]) for k in range(N)), 1.0)
        result = sympy.Integer(0)
        for nn in itertools.product(*[range(pk + 1) for pk in p]):
            pn = tuple(p[k] - nn[k] for k in range(N))
            split_moment = functools.reduce(operator.mul,
                                            (self._in_vars[k].moment(nn[k]) *
                                             self._in_vars[k].moment(pn[k]) for k in range(N)), 1.0)
            result = result + self._coeffs[nn] * self._coeffs[pn] * (full_moment - split_moment)
        return dev_prod * result

    def varOrder(self, n: int) -> sympy.Expr:
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise TaylorException(f'n must be a positive int, got {n}')
        if n > self._max_order:
            raise TaylorException(f'n {n} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        return sum((self.varAt(*p) for p in _multi_indices_with_min1(N, n)), sympy.Integer(0))

    def biasAt(self, *orders: int) -> sympy.Expr:
        if len(orders) != len(self._in_vars):
            raise TaylorException(
                f'expected {len(self._in_vars)} orders, got {len(orders)}')
        for o in orders:
            if not isinstance(o, int) or isinstance(o, bool) or o < 0:
                raise TaylorException(f'each order must be a non-negative int, got {o}')
        total = sum(orders)
        if total < 1:
            raise TaylorException(f'at least one order must be positive, got {orders}')
        if total > self._max_order:
            raise TaylorException(f'total order {total} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        p = orders
        dev_prod = functools.reduce(operator.mul,
                                    (self._in_vars[k].deviation ** p[k] for k in range(N)),
                                    sympy.Integer(1))
        moment_prod = functools.reduce(operator.mul,
                                       (self._in_vars[k].moment(p[k]) for k in range(N)), 1.0)
        return dev_prod * self._coeffs[p] * moment_prod

    def biasOrder(self, n: int) -> sympy.Expr:
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise TaylorException(f'n must be a positive int, got {n}')
        if n > self._max_order:
            raise TaylorException(f'n {n} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        return sum((self.biasAt(*p) for p in _multi_indices(N, n)), sympy.Integer(0))

