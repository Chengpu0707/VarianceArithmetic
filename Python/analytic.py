"""Symbolic statistical Taylor expansion: bounded moments ζ(n, κ), input random
variables (InVar), and the Taylor coefficient / mean (bias) / variance machinery.

Expression conventions:
- coeffs[(α₁,…,α_N)] = (1/α!) · ∂^|α| f / ∏ ∂x_k^{α_k}  (i.e. includes the 1/n! factor).
- moment(n) returns the bound moment ζ(n, κ): symbolic for Gaussian, numeric for Uniform.
- All distributions are assumed symmetric, so ζ(odd, κ) = 0.
"""

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
        # Symmetric distribution assumption: odd-order moments vanish identically.
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
                 distr_type: EDistrType,
                 kappa: typing.Union[float, sympy.Symbol] = None,
                 samples: int = 10000):
        """An input random variable with `value` ± `deviation`, drawn from
        `distr_type` and bounded by ±`kappa`·deviation. Defaults: kappa=5.0
        for Gaussian (leakage ~5e-7), kappa=√3 for Uniform (the maximum
        bounding range of the Uniform distribution). For Uniform, kappa may
        also be a sympy.Symbol — the moment formulas then evaluate symbolically.
        Raises InVarException if any argument has the wrong type or is out of range."""
        if not isinstance(value, sympy.Symbol):
            raise InVarException(f'value must be a sympy.Symbol, got {type(value)}')
        if not isinstance(deviation, sympy.Symbol):
            raise InVarException(f'deviation must be a sympy.Symbol, got {type(deviation)}')
        if not isinstance(distr_type, EDistrType):
            raise InVarException(f'distr_type must be an EDistrType, got {type(distr_type)}')
        if kappa is None:
            # Gaussian default κ=5 keeps bounding leakage ~5e-7.
            # Uniform default κ=√3 makes ζ(0)=ζ(2)=1 exactly (the standard normalization).
            kappa = 5.0 if distr_type == EDistrType.Gaussian else math.sqrt(3)
        if isinstance(kappa, sympy.Symbol):
            # Symbolic κ: only Uniform supports it. Gaussian needs a numeric κ
            # because zeta(n, κ) is evaluated against the precomputed moment table.
            if distr_type != EDistrType.Uniform:
                raise InVarException(
                    f'symbolic kappa is only supported for Uniform, got {distr_type}')
            # Reject if sympy assumptions imply kappa is outside (0, √3].
            if kappa.is_positive is False:
                raise InVarException(f'kappa must be positive, got {kappa}')
            if (kappa - sympy.sqrt(3)).is_positive:
                raise InVarException(
                    f'kappa must be <= sqrt(3) for Uniform distribution, got {kappa}')
        elif isinstance(kappa, float):
            if kappa <= 0:
                raise InVarException(f'kappa must be positive, got {kappa}')
            if distr_type == EDistrType.Uniform and kappa > math.sqrt(3):
                raise InVarException(f'kappa must be <= sqrt(3) for Uniform distribution, got {kappa}')
        else:
            raise InVarException(
                f'kappa must be a float or (for Uniform) sympy.Symbol, got {type(kappa)}')
        if not isinstance(samples, int) or samples <= 0:
            raise InVarException(f'samples must be a positive int, got {samples}')
        self._value = value
        self._deviation = deviation
        self._distr_type = distr_type
        self._kappa = kappa
        self._samples = samples

    @property
    def value(self) -> sympy.Symbol:
        """The sympy symbol for the variable's mean (e.g. `x`)."""
        return self._value

    @property
    def deviation(self) -> sympy.Symbol:
        """The sympy symbol for the variable's standard deviation (e.g. `dx`)."""
        return self._deviation

    @property
    def distr_type(self) -> EDistrType:
        """The underlying distribution: Gaussian or Uniform."""
        return self._distr_type

    @property
    def kappa(self) -> typing.Union[float, sympy.Symbol]:
        """The bounding factor: the variable is restricted to mean ± κ·deviation.
        A float for Gaussian; either a float or a sympy.Symbol for Uniform."""
        return self._kappa

    @property
    def samples(self) -> int:
        """Sample-count metadata (used for downstream variance estimation)."""
        return self._samples

    def moment(self, order: int) -> typing.Union[float, typing.Callable[[int, float], float]]:
        """Return the bound moment ζ(n, κ) of the normalized distribution.
        For Gaussian: a symbolic `zeta(order, kappa)` (numeric via `.evalf()`).
        For Uniform: a Python float per Formula (2.22), or 0 for odd order."""
        if self._distr_type == EDistrType.Gaussian:
            return zeta(order, self._kappa)
        if order % 2 == 1:
            return 0.0
        # Formula (2.22): ζ(n) = 2ρ(κ) · κ^(n+1) / (n+1), 2ρ(κ) = 1/√3
        # Use sympy.sqrt(3) (not math.sqrt(3)) so the result stays exact symbolic
        # for symbolic kappa, and Float-coefficient FP residuals don't appear.
        return self._kappa ** (order + 1) / (sympy.sqrt(3) * (order + 1))


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
        """Build the Taylor coefficient table for `function` around the InVars in
        `in_vars`, up to total order `max_order`. Raises TaylorException if the
        function references any deviation symbol from `in_vars`."""
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
        # Recurrence c[α] = diff(c[α - e_i], x_i) / α_i for any i with α_i > 0.
        # Equivalent to (1/α!) · ∂^|α| f / ∏ ∂x_k^{α_k} but reuses lower-order results.
        for k in range(1, max_order + 1):
            for alpha in _multi_indices(n_vars, k):
                for i, ai in enumerate(alpha):
                    if ai > 0:
                        parent = alpha[:i] + (ai - 1,) + alpha[i + 1:]
                        self._coeffs[alpha] = sympy.diff(self._coeffs[parent], symbols[i]) / ai
                        break

    @property
    def function(self) -> sympy.Expr:
        """The original symbolic function f passed at construction."""
        return self._function

    @property
    def in_vars(self) -> tuple:
        """Tuple of InVar inputs, one per independent variable."""
        return self._in_vars

    @property
    def max_order(self) -> int:
        """Largest total derivative order |α| stored in `coeffs`."""
        return self._max_order

    @property
    def coeffs(self) -> dict:
        """Dict {α: (1/α!) · ∂^|α| f / ∏ ∂x_k^{α_k}} for every multi-index α with
        |α| ≤ max_order. Keys are tuples of length len(in_vars)."""
        return self._coeffs

    def at(self, *orders: int) -> sympy.Expr:
        """Return the Taylor coefficient at multi-index `orders` — i.e.
        coeffs[orders] = (1/α!) · ∂^|α| f / ∏ ∂x_k^{α_k}. Each order must be
        a non-negative int and sum(orders) ≤ max_order."""
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
        """Variance contribution at multi-index p=orders: one term of δ²f as a
        sympy expression. Each p_k must be a non-negative int, sum(p) ≥ 1, and
        sum(p) ≤ max_order. Implements one term of the variance per Formula
        (2.7) (1D) / (2.10) (2D) / extension to N-D. Allows any p_k = 0 as long
        as |p| ≥ 1, which captures 1D contributions (e.g. δ²(x+y) gets nonzero
        terms at p=(2,0) and (0,2))."""
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
        # The j=0 and j=|p| terms (where one factor is f itself) contribute
        # ζ(p)·(1 - ζ(0)), which vanishes when ζ(0)=1 (e.g. Uniform default).
        dev_prod = functools.reduce(operator.mul,
                                    (self._in_vars[k].deviation ** p[k] for k in range(N)),
                                    sympy.Integer(1))
        full_moment = functools.reduce(operator.mul,
                                       (self._in_vars[k].moment(p[k]) for k in range(N)),
                                       sympy.Integer(1))
        result = sympy.Integer(0)
        for nn in itertools.product(*[range(pk + 1) for pk in p]):
            pn = tuple(p[k] - nn[k] for k in range(N))
            split_moment = functools.reduce(operator.mul,
                                            (self._in_vars[k].moment(nn[k]) *
                                             self._in_vars[k].moment(pn[k]) for k in range(N)),
                                            sympy.Integer(1))
            result = result + self._coeffs[nn] * self._coeffs[pn] * (full_moment - split_moment)
        return dev_prod * result

    def varOrder(self, n: int) -> sympy.Expr:
        """Sum of varAt(*p) over all multi-indices p with |p| = n. The total
        variance δ²f equals sum(varOrder(n) for n in 1..max_order)."""
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise TaylorException(f'n must be a positive int, got {n}')
        if n > self._max_order:
            raise TaylorException(f'n {n} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        return sum((self.varAt(*p) for p in _multi_indices(N, n)), sympy.Integer(0))

    def biasAt(self, *orders: int) -> sympy.Expr:
        """Bias (mean − f) contribution at multi-index p=orders: one term of the
        Taylor expansion of the mean as a sympy expression. Each p_k must be a
        non-negative int, sum(p) ≥ 1, and sum(p) ≤ max_order. Implements one
        term of the mean per Formula (2.6) (1D) / (2.9) (2D) / extension to N-D.
        The all-zero p is excluded because that term equals f itself (modulo
        ζ(0)) rather than a bias contribution."""
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
                                       (self._in_vars[k].moment(p[k]) for k in range(N)),
                                       sympy.Integer(1))
        return dev_prod * self._coeffs[p] * moment_prod

    def biasOrder(self, n: int) -> sympy.Expr:
        """Sum of biasAt(*p) over all multi-indices p with |p| = n. The total
        bias (mean − f) equals sum(biasOrder(n) for n in 1..max_order)."""
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise TaylorException(f'n must be a positive int, got {n}')
        if n > self._max_order:
            raise TaylorException(f'n {n} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        return sum((self.biasAt(*p) for p in _multi_indices(N, n)), sympy.Integer(0))

