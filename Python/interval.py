"""Best-effort double-precision interval arithmetic with ULP outward padding
after each operation. Not IEEE-1788 rigorous (no directed-rounding mode),
but tight enough for diagnostic comparison against variance arithmetic.

Conversion from VarDbl: [value - max(unc, ulp(value)), value + max(unc, ulp(value))],
so degenerate (unc=0) sources like IndexSin Quart still get a proper ULP-wide interval.
"""
import math
import sys
import typing

import varDbl


def _safe_ulp(x: float) -> float:
    u = math.ulp(x)
    return u if u > 0.0 else sys.float_info.min


def _padded(lo: float, hi: float) -> 'Interval':
    return Interval(lo - _safe_ulp(lo), hi + _safe_ulp(hi))


class Interval:
    __slots__ = ('_lo', '_hi')

    def __init__(self, lo: float, hi: typing.Optional[float] = None):
        if hi is None:
            hi = lo
        if lo > hi:
            lo, hi = hi, lo
        self._lo = float(lo)
        self._hi = float(hi)

    @classmethod
    def from_varDbl(cls, v: varDbl.VarDbl) -> 'Interval':
        """Half-width = max(uncertainty, ulp(value))."""
        w = max(v.uncertainty(), _safe_ulp(v.value()))
        return cls(v.value() - w, v.value() + w)

    @classmethod
    def centered(cls, value: float, unc: float) -> 'Interval':
        if unc < 0:
            raise ValueError('Interval.centered: unc must be >= 0')
        w = max(unc, _safe_ulp(value))
        return cls(value - w, value + w)

    def lo(self) -> float:   return self._lo
    def hi(self) -> float:   return self._hi
    def mid(self) -> float:  return 0.5 * (self._lo + self._hi)
    def rad(self) -> float:  return 0.5 * (self._hi - self._lo)

    def contains(self, x: float) -> bool:
        return self._lo <= x <= self._hi

    def encloses(self, other: 'Interval') -> bool:
        return self._lo <= other._lo and other._hi <= self._hi

    def __neg__(self) -> 'Interval':
        return Interval(-self._hi, -self._lo)

    def __add__(self, other):
        if isinstance(other, Interval):
            return _padded(self._lo + other._lo, self._hi + other._hi)
        return _padded(self._lo + other, self._hi + other)

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, Interval):
            return _padded(self._lo - other._hi, self._hi - other._lo)
        return _padded(self._lo - other, self._hi - other)

    def __rsub__(self, other):
        return _padded(other - self._hi, other - self._lo)

    def __mul__(self, other):
        if isinstance(other, Interval):
            a = self._lo * other._lo; b = self._lo * other._hi
            c = self._hi * other._lo; d = self._hi * other._hi
            return _padded(min(a, b, c, d), max(a, b, c, d))
        a = self._lo * other; b = self._hi * other
        return _padded(min(a, b), max(a, b))

    __rmul__ = __mul__

    def __repr__(self) -> str:
        return f'[{self._lo!r}, {self._hi!r}]'

    def __str__(self) -> str:
        return f'[{self._lo}, {self._hi}]'
