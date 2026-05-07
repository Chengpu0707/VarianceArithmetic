"""Symbolic statistical Taylor expansion: bounded moments ζ(n, κ), input random
variables (InVar), and the Taylor coefficient / mean (bias) / variance machinery.

Expression conventions:
- coeffs[(α₁,…,α_N)] = (1/α!) · ∂^|α| f / ∏ ∂x_k^{α_k}  (i.e. includes the 1/n! factor).
- moment(n) returns the bound moment ζ(n, κ): symbolic for Gaussian, numeric for Uniform.
- All distributions are assumed symmetric, so ζ(odd, κ) = 0.
"""

import csv
import enum
import functools
import itertools
import math
import operator
import sympy
import typing

import moment


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
            return sympy.Float(moment.Normal(bounding=float(kappa))[int(n)])


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
            return 0
        # Formula (2.2) normalized: ζ(n, κ) = ∫z^n ρ dz / ∫ρ dz.
        # For Uniform on [-√3, √3] truncated to [-κ, κ]: ζ(n, κ) = κ^n / (n+1).
        # ζ(0, κ) = 1 by construction (normalized).
        return self._kappa ** order / (order + 1)


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


class StatTaylor:
    __slots__ = ('_function', '_in_vars', '_max_order', '_coeffs')

    def __init__(self, function: sympy.Expr, in_vars: tuple, max_order: int = 16):
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
        # Optimization: if the function is identically zero, every higher-order
        # derivative is also zero — skip the (potentially huge) recurrence
        # loop. This lets StatMatrix use any max_order for its placeholder
        # Integer(0) function without building a multi-million-entry coeff
        # table; lookups in at()/varAt()/biasAt() use .get(..., 0) so missing
        # entries report as zero.
        if function == 0:
            return
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
        return self._coeffs.get(orders, sympy.Integer(0))

    def varAt(self, *orders: int) -> sympy.Expr:
        """Variance contribution at multi-index p=orders: one term of δ²f as a
        sympy expression per Formula (2.7) (1D) / (2.10) (2D) / extension to N-D.
        Each p_k must be a non-negative int with sum(p) ≤ max_order. The all-zero
        p=(0,…,0) is now valid: it returns f(x)²·∏ζ_k(0)·(1−∏ζ_k(0)), the n=0
        contribution that vanishes only when ζ(0)=1."""
        if len(orders) != len(self._in_vars):
            raise TaylorException(
                f'expected {len(self._in_vars)} orders, got {len(orders)}')
        for o in orders:
            if not isinstance(o, int) or isinstance(o, bool) or o < 0:
                raise TaylorException(f'each order must be a non-negative int, got {o}')
        total = sum(orders)
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
            result = result + (self._coeffs.get(nn, sympy.Integer(0))
                               * self._coeffs.get(pn, sympy.Integer(0))
                               * (full_moment - split_moment))
        return dev_prod * result

    def varOrder(self, n: int) -> sympy.Expr:
        """Sum of varAt(*p) over all multi-indices p with |p| = n. The total
        variance δ²f equals sum(varOrder(n) for n in 0..max_order)."""
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise TaylorException(f'n must be a non-negative int, got {n}')
        if n > self._max_order:
            raise TaylorException(f'n {n} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        return sum((self.varAt(*p) for p in _multi_indices(N, n)), sympy.Integer(0))

    def biasAt(self, *orders: int) -> sympy.Expr:
        """One term of the mean's Taylor expansion at multi-index p=orders, per
        Formula (2.6) (1D) / (2.9) (2D) / extension to N-D. Each p_k must be a
        non-negative int with sum(p) ≤ max_order. The all-zero p=(0,…,0) returns
        f(x)·∏ζ_k(0); the bias `mean − f` therefore equals (sum of biasAt over
        all p) − f, with the (0,…,0) term contributing f·(∏ζ_k(0) − 1)."""
        if len(orders) != len(self._in_vars):
            raise TaylorException(
                f'expected {len(self._in_vars)} orders, got {len(orders)}')
        for o in orders:
            if not isinstance(o, int) or isinstance(o, bool) or o < 0:
                raise TaylorException(f'each order must be a non-negative int, got {o}')
        total = sum(orders)
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
        return dev_prod * self._coeffs.get(p, sympy.Integer(0)) * moment_prod

    def biasOrder(self, n: int) -> sympy.Expr:
        """Sum of biasAt(*p) over all multi-indices p with |p| = n. The total
        mean equals sum(biasOrder(n) for n in 0..max_order); the bias is the
        total mean minus f."""
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise TaylorException(f'n must be a non-negative int, got {n}')
        if n > self._max_order:
            raise TaylorException(f'n {n} exceeds max_order {self._max_order}')
        N = len(self._in_vars)
        return sum((self.biasAt(*p) for p in _multi_indices(N, n)), sympy.Integer(0))

    def dump(self, path: str) -> None:
        """Write varAt(*p) and biasAt(*p) for every multi-index p with
        sum(p) ≤ max_order to `path` in CSV format. Rows where both varAt
        and biasAt are 0 are skipped.

        The first row is a single-field `function: <expr>` line with each
        InVar's value symbol `x` rendered as `x~dx` (where `dx` is its
        deviation symbol). The second row is the column header:
        `order` + one column per InVar (named `x~dx`) + `varAt` + `biasAt`.
        Each subsequent row holds one multi-index (the per-InVar columns
        give the derivative order per input variable), in lexicographic
        order within each total-order block. Field quoting is handled by
        the standard csv writer (commas in the function expression are
        quoted automatically)."""
        if not isinstance(path, str):
            raise TaylorException(f'path must be a str, got {type(path)}')
        N = len(self._in_vars)
        # Display each InVar value symbol `x` as `x~dx` in the function and
        # use the same `x~dx` strings as per-InVar column headers.
        invar_labels = [f'{inv.value}~{inv.deviation}' for inv in self._in_vars]
        fn_subs = {inv.value: sympy.Symbol(label)
                   for inv, label in zip(self._in_vars, invar_labels)}
        fn_display = self._function.subs(fn_subs)
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([f'function: {fn_display}'])
            writer.writerow(['order'] + invar_labels + ['varAt', 'biasAt'])
            for n in range(self._max_order + 1):
                for p in _multi_indices(N, n):
                    v = self.varAt(*p)
                    b = self.biasAt(*p)
                    if v == 0 and b == 0:
                        continue
                    writer.writerow([n] + list(p) + [str(v), str(b)])


class StatMatrix(StatTaylor):
    """N×N matrix whose entries are a mix of InVar (uncertain) and value
    (certain) cells. Inherits StatTaylor with a placeholder zero function
    (max_order=0); the symbolic determinant and per-entry StatTaylor are
    exposed via determ() and item(). Indices are row-major: index(row, col)
    = row*N + col."""

    __slots__ = ('_N', '_matrix', '_invar_pos')

    def __init__(self, N: int, items: dict,
                 in_vars: typing.Optional[tuple] = None,
                 max_order: int = 16):
        """`N`: positive int, side length.
        `items`: dict {(row, col): InVar | numeric/sympy.Expr}. Any (row, col)
        not present in `items` defaults to symbolic zero — `sympy.Integer(0)`,
        which is the `sympy.S.Zero` singleton — so the resulting `matrix` cell
        is always a sympy expression.
        `in_vars`: optional explicit tuple of InVar to use as the underlying
        StatTaylor inputs. If None (default), inferred from InVar entries in
        `items` (and at least one such entry is required). If provided (as for
        derived matrices like adjugate(), where every entry is a sympy
        expression in the original in_vars rather than an InVar), used directly;
        any InVar entries in `items` are then stored as their value symbol.
        `max_order`: forwarded to the underlying StatTaylor (default 16) and
        used as the default for derived analyses (item, dump). The parent
        StatTaylor's coefficient table is never materialised because the
        StatMatrix's function is the placeholder Integer(0); see the
        zero-function shortcut in StatTaylor.__init__."""
        if not isinstance(N, int) or isinstance(N, bool) or N <= 0:
            raise TaylorException(f'N must be a positive int, got {N}')
        if not isinstance(items, dict):
            raise TaylorException(f'items must be a dict, got {type(items)}')
        for key in items:
            if not (isinstance(key, tuple) and len(key) == 2):
                raise TaylorException(
                    f'items key must be a (row, col) tuple, got {key}')
            r, c = key
            if (not isinstance(r, int) or isinstance(r, bool)
                    or not isinstance(c, int) or isinstance(c, bool)):
                raise TaylorException(f'items key {key} must contain ints')
            if not (0 <= r < N and 0 <= c < N):
                raise TaylorException(f'items key {key} out of range for N={N}')
        if in_vars is not None:
            if not isinstance(in_vars, tuple):
                raise TaylorException(
                    f'in_vars must be a tuple, got {type(in_vars)}')
            if not in_vars:
                raise TaylorException('in_vars must have at least one element')
            for v in in_vars:
                if not isinstance(v, InVar):
                    raise TaylorException(
                        f'each in_vars element must be InVar, got {type(v)}')
        item_in_vars = []
        invar_pos = {}  # (r, c) -> position in item_in_vars (only used when in_vars is None)
        rows = []
        for r in range(N):
            row = []
            for c in range(N):
                entry = items.get((r, c), sympy.Integer(0))
                if isinstance(entry, InVar):
                    if in_vars is None:
                        invar_pos[(r, c)] = len(item_in_vars)
                        item_in_vars.append(entry)
                    row.append(entry.value)
                else:
                    row.append(sympy.sympify(entry))
            rows.append(row)
        if in_vars is None:
            if not item_in_vars:
                raise TaylorException(
                    'StatMatrix must contain at least one InVar entry '
                    '(or pass in_vars explicitly for a derived matrix)')
            in_vars = tuple(item_in_vars)
        self._N = N
        self._matrix = sympy.Matrix(rows)
        self._invar_pos = invar_pos
        super().__init__(sympy.Integer(0), in_vars, max_order=max_order)

    @property
    def N(self) -> int:
        """Side length of the square matrix."""
        return self._N

    @property
    def matrix(self) -> sympy.Matrix:
        """sympy.Matrix with value-symbols at InVar positions and sympy
        expressions / numerics at all other positions."""
        return self._matrix

    def index(self, row: int, col: int) -> int:
        """Row-major flat index row*N + col for the entry at (row, col).
        For all-InVar matrices (e.g. WorstMatrix) this is also the index
        into `in_vars`; for mixed matrices it is just the layout position."""
        if (not isinstance(row, int) or isinstance(row, bool)
                or not isinstance(col, int) or isinstance(col, bool)):
            raise TaylorException(f'row and col must be ints, got ({row}, {col})')
        if not (0 <= row < self._N and 0 <= col < self._N):
            raise TaylorException(
                f'(row, col)=({row}, {col}) out of range for N={self._N}')
        return row * self._N + col

    def pos(self, index: int) -> tuple:
        """Inverse of index(): (row, col) for the given flat index."""
        if not isinstance(index, int) or isinstance(index, bool):
            raise TaylorException(f'index must be an int, got {index}')
        if not (0 <= index < self._N * self._N):
            raise TaylorException(
                f'index={index} out of range for N={self._N}')
        return (index // self._N, index % self._N)

    def subMatrix(self, positions) -> 'StatMatrix':
        """Return a new StatMatrix with every row and every col appearing in
        `positions` removed. Each element of `positions` is a (row, col) tuple.
        Unique rows and unique cols must be equinumerous so the result stays
        square; surviving InVar instances are reused (same symbols), and
        non-zero value entries are preserved."""
        if not isinstance(positions, (list, tuple)):
            raise TaylorException(
                f'positions must be a list/tuple, got {type(positions)}')
        rows_to_drop = set()
        cols_to_drop = set()
        for p in positions:
            if not (isinstance(p, tuple) and len(p) == 2):
                raise TaylorException(f'each position must be a 2-tuple, got {p}')
            r, c = p
            if (not isinstance(r, int) or isinstance(r, bool)
                    or not isinstance(c, int) or isinstance(c, bool)):
                raise TaylorException(f'position must be ints, got ({r}, {c})')
            if not (0 <= r < self._N and 0 <= c < self._N):
                raise TaylorException(
                    f'position ({r}, {c}) out of range for N={self._N}')
            rows_to_drop.add(r)
            cols_to_drop.add(c)
        if len(rows_to_drop) != len(cols_to_drop):
            raise TaylorException(
                f'subMatrix needs equal row/col deletions; got '
                f'{len(rows_to_drop)} unique rows and {len(cols_to_drop)} unique cols')
        new_N = self._N - len(rows_to_drop)
        if new_N <= 0:
            raise TaylorException('subMatrix would empty the matrix')
        survivor_rows = [r for r in range(self._N) if r not in rows_to_drop]
        survivor_cols = [c for c in range(self._N) if c not in cols_to_drop]
        new_items = {}
        for new_r, r in enumerate(survivor_rows):
            for new_c, c in enumerate(survivor_cols):
                if (r, c) in self._invar_pos:
                    new_items[(new_r, new_c)] = self._in_vars[self._invar_pos[(r, c)]]
                else:
                    val = self._matrix[r, c]
                    if val != 0:
                        new_items[(new_r, new_c)] = val
        # If no InVar items survive (e.g. subMatrix of an adjugate), inherit
        # the original's in_vars so the derived StatMatrix stays valid.
        if any(isinstance(v, InVar) for v in new_items.values()):
            return StatMatrix(new_N, new_items, max_order=self._max_order)
        return StatMatrix(new_N, new_items, in_vars=self._in_vars,
                          max_order=self._max_order)

    def determ(self, max_order: int = None) -> StatTaylor:
        """Return a fresh StatTaylor for the symbolic determinant of this matrix.
        The default max_order is 2*N — sufficient for variance analysis of the
        polynomial determinant. Reuses the matrix's in_vars."""
        if max_order is None:
            max_order = 2 * self._N
        det_expr = self._matrix.det()
        return StatTaylor(det_expr, self._in_vars, max_order=max_order)

    def adjugate(self) -> 'StatMatrix':
        """Return the adjugate (classical adjoint) as a new StatMatrix. Each
        entry is the corresponding cofactor expression (sympy.Expr in the
        original in_vars):
            adj(M)[i, j] = (-1)^(i+j) · det(M with row j and col i removed)
        Equivalently, adj(M) is the transpose of the cofactor matrix. The
        defining identity is `adj(M) · M == det(M) · I`. The returned matrix
        has no InVar items of its own; it inherits this matrix's `in_vars`."""
        adj_matrix = self._matrix.adjugate()
        items = {}
        for i in range(self._N):
            for j in range(self._N):
                entry = adj_matrix[i, j]
                if entry != 0:
                    items[(i, j)] = entry
        return StatMatrix(self._N, items, in_vars=self._in_vars,
                          max_order=self._max_order)

    def reverse(self) -> 'StatMatrix':
        """Return the matrix inverse as a new StatMatrix. Each entry is the
        rational expression M^(-1)[i, j] = adj(M)[i, j] / det(M). The defining
        identity is `reverse(M) · M == I == M · reverse(M)`. The returned
        matrix has no InVar items of its own; it inherits this matrix's
        `in_vars`. Raises sympy's NonInvertibleMatrixError if the matrix is
        singular (det = 0)."""
        inv_matrix = self._matrix.inv()
        items = {}
        for i in range(self._N):
            for j in range(self._N):
                entry = inv_matrix[i, j]
                if entry != 0:
                    items[(i, j)] = entry
        return StatMatrix(self._N, items, in_vars=self._in_vars,
                          max_order=self._max_order)

    def dump(self, path: str, max_order: int = None) -> None:
        """Write per-entry varAt and biasAt to `path` in CSV format
        (overrides StatTaylor.dump for matrices). For each (row, col) of
        the matrix, builds a StatTaylor of the entry over all in_vars and
        emits one row per nonzero (multi-index, varAt, biasAt) triple,
        tagged with (row, col) so all entries share one file. Defaults to
        this StatMatrix's own max_order (16 unless overridden at construction).

        File structure (CSV):
          - One single-field `function (row, col): <expr>` row per entry,
            with each InVar's value symbol rendered as `x~dx` (deviation
            appended). The csv writer quotes the field automatically since
            it contains commas.
          - Header row: `order, row, col` + one column per InVar (named
            `x~dx`) + `varAt`, `biasAt`.
          - Data rows: `order` is the multi-index sum; `row` and `col`
            identify the matrix entry; the per-InVar columns hold the
            multi-index components; the last two columns are varAt and biasAt.
        Rows where both varAt and biasAt are 0 are skipped."""
        if not isinstance(path, str):
            raise TaylorException(f'path must be a str, got {type(path)}')
        if max_order is None:
            max_order = self._max_order
        if (not isinstance(max_order, int) or isinstance(max_order, bool)
                or max_order < 0):
            raise TaylorException(
                f'max_order must be a non-negative int, got {max_order}')
        N = self._N
        n_in_vars = len(self._in_vars)
        invar_labels = [f'{inv.value}~{inv.deviation}' for inv in self._in_vars]
        fn_subs = {inv.value: sympy.Symbol(label)
                   for inv, label in zip(self._in_vars, invar_labels)}
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            for r in range(N):
                for c in range(N):
                    display = self._matrix[r, c].subs(fn_subs)
                    writer.writerow([f'function ({r}, {c}): {display}'])
            writer.writerow(['order', 'row', 'col'] + invar_labels +
                            ['varAt', 'biasAt'])
            for r in range(N):
                for c in range(N):
                    T = StatTaylor(self._matrix[r, c], self._in_vars,
                                   max_order=max_order)
                    for n in range(max_order + 1):
                        for p in _multi_indices(n_in_vars, n):
                            v = T.varAt(*p)
                            b = T.biasAt(*p)
                            if v == 0 and b == 0:
                                continue
                            writer.writerow([n, r, c] + list(p)
                                            + [str(v), str(b)])

    def item(self, position, max_order: int = None) -> StatTaylor:
        """Return a StatTaylor for the entry at `position`=(row, col).
        If the entry is an InVar, returns a single-variable StatTaylor over
        that InVar (function = its value symbol). Otherwise the entry is a
        sympy expression in the original in_vars (e.g. for derived matrices
        produced by adjugate() or reverse()), and the returned StatTaylor
        uses all of this StatMatrix's in_vars. Defaults to this StatMatrix's
        own max_order (16 unless overridden at construction)."""
        if max_order is None:
            max_order = self._max_order
        if not (isinstance(position, tuple) and len(position) == 2):
            raise TaylorException(
                f'position must be a (row, col) tuple, got {position}')
        row, col = position
        # Validate (row, col) range via index().
        self.index(row, col)
        if (row, col) in self._invar_pos:
            inv = self._in_vars[self._invar_pos[(row, col)]]
            return StatTaylor(inv.value, (inv,), max_order=max_order)
        return StatTaylor(self._matrix[row, col], self._in_vars,
                          max_order=max_order)


class WorstMatrix(StatMatrix):
    """All-InVar special case of StatMatrix: every entry is an InVar with value
    symbol `m_{r}_{c}` and deviation symbol `dm_{r}_{c}`, drawn from
    `distr_type` (default Uniform). Used for worst-case variance analysis."""

    __slots__ = ()

    def __init__(self, N: int,
                 distr_type: EDistrType = EDistrType.Uniform,
                 kappa: typing.Union[float, sympy.Symbol] = None,
                 samples: int = 10000,
                 max_order: int = 16):
        if not isinstance(N, int) or isinstance(N, bool) or N <= 0:
            raise TaylorException(f'N must be a positive int, got {N}')
        items = {}
        for r in range(N):
            for c in range(N):
                v = sympy.Symbol(f'm_{r}_{c}')
                d = sympy.Symbol(f'dm_{r}_{c}')
                items[(r, c)] = InVar(v, d, distr_type, kappa=kappa, samples=samples)
        super().__init__(N, items, max_order=max_order)

