"""Unit tests for analytic.py (statistical Taylor expansion).

Validation tests:
  TestInVar              — InVar construction, defaults, immutability, type checks.
  TestInVarMoment        — InVar.moment() for Gaussian (symbolic ζ) and Uniform (numeric).
  TestStatTaylor         — StatTaylor construction, properties, immutability, type checks.
  TestStatTaylorMethod   — StatTaylor.at() argument validation and small structural cases.
  TestStatTaylorVarAt    — varAt(*orders) argument validation and 1D/2D structural cases.
  TestStatTaylorVarOrder — varOrder(n) argument validation and structural sums.
  TestStatTaylorBiasAt   — biasAt(*orders) argument validation (incl. all-zero rejection).
  TestStatTaylorBiasOrder — biasOrder(n) argument validation.

Function-specific tests (all use Uniform default κ=√3 so ζ(0)=1, ζ(odd)=0).
The 1D classes also verify Formulas (2.14)–(2.21) of the Short paper:
  TestLinear             — f = x:               trivial coeffs/var/bias.
  TestQuadratic          — f = x²:              all bias/var collapse to single terms.
  TestExp                — f = e^x:             (2.14) exp mean, (2.15) exp precision.
  TestLog                — f = log(x):          (2.16) log mean, (2.17) log precision.
  TestSine               — f = sin(x):          (2.18) sin mean, (2.19) sin precision.
  TestPow                — f = x^c (c=1.5):     (2.20) power mean, (2.21) power precision.
  TestSinXdivX           — f = sin(x)/x:        rational mixture of sin/cos/powers.

Multi-variable tests (Uniform on every variable):
  TestXaddY              — f = x + y:           Formula (2.10) bias, (2.11) variance.
  TestXmulY              — f = x·y:             Formula (2.12) bias, (2.13) variance.
  TestXdivY              — f = x/y:             linear in x; rational/alternating in y.
  TestXpowY              — f = x^y:             non-separable, mixes x^y and log(x).
  TestXaddYaddZ          — f = x + y + z:       3D extension of TestXaddY.
"""

import csv
import math
import os
import sympy
import unittest

import analytic
from indexSin import OUTDIR


class TestInVar(unittest.TestCase):

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.k = sympy.Symbol('k')

    def test_gaussian_defaults(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.assertEqual(v.value, self.x)
        self.assertEqual(v.deviation, self.dx)
        self.assertEqual(v.distr_type, analytic.EDistrType.Gaussian)
        self.assertEqual(v.kappa, 5.0)
        self.assertEqual(v.samples, 10000)

    def test_uniform_defaults(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertEqual(v.kappa, math.sqrt(3))
        self.assertEqual(v.samples, 10000)

    def test_explicit_kappa(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=3.0)
        self.assertEqual(v.kappa, 3.0)

    def test_explicit_samples(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, samples=500)
        self.assertEqual(v.samples, 500)

    def test_readonly_value(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        with self.assertRaises(AttributeError):
            v.value = sympy.Symbol('y')

    def test_readonly_deviation(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        with self.assertRaises(AttributeError):
            v.deviation = sympy.Symbol('dy')

    def test_readonly_distr_type(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        with self.assertRaises(AttributeError):
            v.distr_type = analytic.EDistrType.Uniform

    def test_readonly_kappa(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        with self.assertRaises(AttributeError):
            v.kappa = 3.0

    def test_readonly_samples(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        with self.assertRaises(AttributeError):
            v.samples = 1

    def test_no_extra_attributes(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        with self.assertRaises(AttributeError):
            v.extra = 42

    def test_invalid_value(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(1.0, self.dx, analytic.EDistrType.Gaussian)

    def test_invalid_deviation(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, 1.0, analytic.EDistrType.Gaussian)

    def test_invalid_distr_type(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, 'Gaussian')

    def test_invalid_kappa_type(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=5)

    def test_invalid_kappa_nonpositive(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=0.0)
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=-1.0)

    def test_invalid_kappa_uniform_exceeds_sqrt3(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=2.0)

    def test_uniform_symbolic_kappa(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=self.k)
        self.assertEqual(v.kappa, self.k)
        self.assertEqual(v.distr_type, analytic.EDistrType.Uniform)

    def test_uniform_symbolic_kappa_unspecified_assumptions(self):
        # Symbol without assumptions: cannot prove out-of-range, accept.
        k = sympy.Symbol('k')
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=k)
        self.assertEqual(v.kappa, k)

    def test_invalid_uniform_symbolic_kappa_negative(self):
        # Symbol with negative=True is provably non-positive → reject.
        k = sympy.Symbol('k', negative=True)
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=k)

    def test_invalid_uniform_symbolic_kappa_zero(self):
        # Symbol with zero=True: is_positive is False → reject.
        k = sympy.Symbol('k', zero=True)
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=k)

    def test_invalid_gaussian_symbolic_kappa(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=self.k)

    def test_invalid_samples_type(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, samples=1.0)

    def test_invalid_samples_nonpositive(self):
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, samples=0)
        with self.assertRaises(analytic.InVarException):
            analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, samples=-1)


class TestInVarMoment(unittest.TestCase):

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.k = sympy.Symbol('k', positive=True)

    def test_gaussian_default_kappa_order0(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.assertEqual(v.moment(0), analytic.zeta(0, v.kappa))

    def test_gaussian_default_kappa_order2(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.assertEqual(v.moment(2), analytic.zeta(2, v.kappa))

    def test_gaussian_odd_order_is_zero(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.assertEqual(v.moment(1), 0)
        self.assertEqual(v.moment(3), 0)

    def test_uniform_default_kappa_order0(self):
        import moment as _mmt
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertAlmostEqual(v.moment(0), _mmt.Uniform(bounding=math.sqrt(3))[0])

    def test_uniform_default_kappa_order2(self):
        import moment as _mmt
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertAlmostEqual(v.moment(2), _mmt.Uniform(bounding=math.sqrt(3))[2])

    def test_uniform_odd_order_is_zero(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertEqual(v.moment(1), 0)

    def test_custom_kappa_gaussian(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=3.0)
        self.assertEqual(v.moment(2), analytic.zeta(2, v.kappa))

    def test_uniform_symbolic_kappa_order0(self):
        # Normalized (2.2): ζ(0, κ) = 1 by construction.
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=self.k)
        self.assertEqual(sympy.simplify(v.moment(0) - 1), 0)

    def test_uniform_symbolic_kappa_order2(self):
        # Normalized: ζ(2, κ) = κ² / 3.
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=self.k)
        self.assertEqual(sympy.simplify(v.moment(2) - self.k**2 / 3), 0)

    def test_uniform_symbolic_kappa_order_odd_is_zero(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=self.k)
        self.assertEqual(v.moment(1), 0)
        self.assertEqual(v.moment(3), 0)


class TestStatTaylor(unittest.TestCase):

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.y = sympy.Symbol('y')
        self.dy = sympy.Symbol('dy')
        self.vx = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.vy = analytic.InVar(self.y, self.dy, analytic.EDistrType.Uniform)

    def test_single_var(self):
        t = analytic.StatTaylor(sympy.sin(self.x), (self.vx,), max_order=2)
        self.assertEqual(t.function, sympy.sin(self.x))
        self.assertEqual(t.in_vars, (self.vx,))

    def test_two_vars(self):
        t = analytic.StatTaylor(self.x + self.y, (self.vx, self.vy), max_order=2)
        self.assertEqual(len(t.in_vars), 2)

    def test_default_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=200)
        self.assertEqual(t.max_order, 200)

    def test_custom_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=5)
        self.assertEqual(t.max_order, 5)

    def test_readonly_function(self):
        t = analytic.StatTaylor(sympy.sin(self.x), (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.function = sympy.cos(self.x)

    def test_readonly_in_vars(self):
        t = analytic.StatTaylor(sympy.sin(self.x), (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.in_vars = (self.vy,)

    def test_readonly_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=5)
        with self.assertRaises(AttributeError):
            t.max_order = 10

    def test_readonly_coeffs(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.coeffs = (self.x,)

    def test_no_extra_attributes(self):
        t = analytic.StatTaylor(sympy.sin(self.x), (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.extra = 42

    def test_invalid_function_type(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor('sin(x)', (self.vx,))

    def test_invalid_in_vars_not_tuple(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(sympy.sin(self.x), [self.vx])

    def test_invalid_in_vars_empty(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(sympy.sin(self.x), ())

    def test_invalid_in_vars_element(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(sympy.sin(self.x), (self.x,))

    def test_invalid_max_order_float(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(self.x, (self.vx,), max_order=2.0)

    def test_invalid_max_order_bool(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(self.x, (self.vx,), max_order=True)

    def test_invalid_max_order_negative(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(self.x, (self.vx,), max_order=-1)


class TestStatTaylorMethod(unittest.TestCase):

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.y = sympy.Symbol('y')
        self.dy = sympy.Symbol('dy')
        self.vx = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.vy = analytic.InVar(self.y, self.dy, analytic.EDistrType.Uniform)

    def _check(self, result, expected):
        self.assertEqual(sympy.expand(result - expected), 0)

    def test_two_vars_order1(self):
        # f(x,y)=xy: at(1,0)=∂f/∂x=y, at(0,1)=∂f/∂y=x
        t = analytic.StatTaylor(self.x*self.y, (self.vx, self.vy), max_order=2)
        self._check(t.at(1, 0), self.y)
        self._check(t.at(0, 1), self.x)

    def test_two_vars_order2(self):
        # f(x,y)=xy: pure partials ∂²f/∂x²=0, ∂²f/∂y²=0 (mixed term not stored)
        t = analytic.StatTaylor(self.x*self.y, (self.vx, self.vy), max_order=2)
        self._check(t.at(2, 0), sympy.Integer(0))
        self._check(t.at(0, 2), sympy.Integer(0))

    def test_invalid_function_contains_deviation(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatTaylor(self.x + self.dx, (self.vx,), max_order=2)

    def test_taylor_exceeds_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(3)

    def test_invalid_order_float(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(1.0)

    def test_invalid_order_bool(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(True)

    def test_invalid_order_negative(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(-1)


class _Base(unittest.TestCase):
    # Symbolic kappa per variable: ζ_k(0) = k_k/√3 ≠ 1 in general, so test
    # expressions must carry the moment factors explicitly.

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.y = sympy.Symbol('y')
        self.dy = sympy.Symbol('dy')
        self.z = sympy.Symbol('z')
        self.dz = sympy.Symbol('dz')
        self.k_x = sympy.Symbol('k_x', positive=True)
        self.k_y = sympy.Symbol('k_y', positive=True)
        self.k_z = sympy.Symbol('k_z', positive=True)
        self.vx = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform, kappa=self.k_x)
        self.vy = analytic.InVar(self.y, self.dy, analytic.EDistrType.Uniform, kappa=self.k_y)
        self.vz = analytic.InVar(self.z, self.dz, analytic.EDistrType.Uniform, kappa=self.k_z)

    def _check(self, result, expected):
        self.assertEqual(sympy.simplify(result - expected), 0)

    # 1D helpers (require self.t set by subclass): moment(odd)=0 for symmetric.

    def _check_varAt_order1_zero(self):
        self._check(self.t.varAt(1), sympy.Integer(0))

    def _check_biasAt_odd_zero(self, n):
        self._check(self.t.biasAt(n), sympy.Integer(0))

    def _check_orders_match_at(self, n):
        # 1D: varOrder(n) == varAt(n), biasOrder(n) == biasAt(n)
        self._check(self.t.varOrder(n), self.t.varAt(n))
        self._check(self.t.biasOrder(n), self.t.biasAt(n))

class TestStatTaylorVarAt(_Base):

    def test_invalid_wrong_num_orders(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(1, 1)

    def test_invalid_order_float(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(2.0)

    def test_invalid_order_bool(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(True)

    def test_invalid_order_exceeds_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(5)

    def test_2d_product_order2_2(self):
        # f(x,y)=x*y: p=(2,2), nn=(1,1) contributes (mixed partials beyond order 1 are 0)
        t = analytic.StatTaylor(self.x * self.y, (self.vx, self.vy), max_order=4)
        self._check(t.varAt(2, 2),
                    self.vx.moment(2) * self.vy.moment(2) * self.dx**2 * self.dy**2)


class TestStatTaylorVarOrder(_Base):

    def test_invalid_n_float(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(2.0)

    def test_invalid_n_bool(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(True)

    def test_invalid_n_exceeds_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(5)

    def test_2d_product_order4(self):
        # f(x,y)=x*y: only varAt(2,2) contributes among (1,3),(2,2),(3,1)
        t = analytic.StatTaylor(self.x * self.y, (self.vx, self.vy), max_order=4)
        self._check(t.varOrder(4),
                    self.vx.moment(2) * self.vy.moment(2) * self.dx**2 * self.dy**2)


class TestStatTaylorBiasAt(_Base):

    def test_invalid_wrong_num_orders(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(1, 1)

    def test_invalid_order_float(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(2.0)

    def test_invalid_order_bool(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(True)

    def test_invalid_order_negative(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(-1)

    def test_invalid_order_exceeds_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(5)

    def test_2d_partial_zero_allowed(self):
        # biasAt(1,0)=δx · ∂f/∂x · ζ_x(1) · ζ_y(0); for f=x and symmetric x, ζ_x(1)=0
        t = analytic.StatTaylor(self.x, (self.vx, self.vy), max_order=4)
        self._check(t.biasAt(1, 0), sympy.Integer(0))


class TestStatTaylorBiasOrder(_Base):

    def test_invalid_n_float(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(2.0)

    def test_invalid_n_bool(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(True)

    def test_invalid_n_exceeds_max_order(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(5)


class TestStatTaylorDump(_Base):

    def test_dump_1d_writes_function_header_and_rows(self):
        import tempfile
        t = analytic.StatTaylor(sympy.sin(self.x), (self.vx,), max_order=3)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'dump.csv')
            t.dump(path)
            with open(path, newline='') as f:
                rows = list(csv.reader(f))
        # Function row: single-field line with InVar value as `x~dx`.
        self.assertEqual(rows[0], ['function: sin(x~dx)'])
        # Header: order + per-InVar columns (named `x~dx`) + varAt + biasAt.
        self.assertEqual(rows[1], ['order', 'x~dx', 'varAt', 'biasAt'])
        # 1D max_order=3 sin(x): orders 0 and 2 are nonzero (1, 3 skipped).
        self.assertEqual(len(rows), 2 + 2)
        # Row 0: order=0, x~dx=0, varAt=0, biasAt=sin(x)
        self.assertEqual(rows[2][0], '0')
        self.assertEqual(rows[2][1], '0')
        self.assertEqual(rows[2][2], str(t.varAt(0)))
        self.assertEqual(rows[2][3], str(t.biasAt(0)))
        # Row 1: order=2, x~dx=2
        self.assertEqual(rows[3][0], '2')
        self.assertEqual(rows[3][1], '2')
        self.assertEqual(rows[3][2], str(t.varAt(2)))
        self.assertEqual(rows[3][3], str(t.biasAt(2)))

    def test_dump_skips_all_zero_rows(self):
        import tempfile
        # 1D linear x: every order ≥ 2 has varAt=biasAt=0; order 1 has both 0
        # (moment(1)=0), order 0 has biasAt=x (nonzero) and varAt=0.
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'dump.csv')
            t.dump(path)
            with open(path, newline='') as f:
                rows = list(csv.reader(f))
        # Expect: function + header + (order=0, biasAt=x) + (order=2, varAt=dx²·ζ(2))
        self.assertEqual(len(rows), 2 + 2)
        self.assertEqual(rows[2][0], '0')
        self.assertEqual(rows[3][0], '2')

    def test_dump_2d_per_invar_columns(self):
        import tempfile
        t = analytic.StatTaylor(self.x * self.y,
                                (self.vx, self.vy), max_order=2)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'dump.csv')
            t.dump(path)
            with open(path, newline='') as f:
                rows = list(csv.reader(f))
        # Header has 5 columns: order, x~dx, y~dy, varAt, biasAt.
        self.assertEqual(rows[1], ['order', 'x~dx', 'y~dy', 'varAt', 'biasAt'])
        # Nonzero rows for x*y at max_order=2:
        # (0,0): biasAt=x*y; (0,2): varAt=dy²·x²; (2,0): varAt=dx²·y²
        self.assertEqual(len(rows), 2 + 3)
        # Each data row's order column should equal sum of per-InVar columns.
        for row in rows[2:]:
            self.assertEqual(int(row[0]), int(row[1]) + int(row[2]))

    def test_dump_invalid_path(self):
        t = analytic.StatTaylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.dump(123)

    def test_dump_2d_function_uses_invar_notation(self):
        import tempfile
        # Function references both InVars; both must render as `x~dx` form.
        t = analytic.StatTaylor(self.x * self.y + sympy.sin(self.x),
                                (self.vx, self.vy), max_order=1)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'dump.csv')
            t.dump(path)
            with open(path, newline='') as f:
                first_row = next(csv.reader(f))
        # Single-field function row.
        self.assertEqual(len(first_row), 1)
        first_line = first_row[0]
        self.assertIn('x~dx', first_line)
        self.assertIn('y~dy', first_line)
        # Function line shouldn't carry bare `x` or `y` outside the `~` token.
        # (Strip the `x~dx`/`y~dy` substrings, then check no leftover `x`/`y`.)
        stripped = first_line.replace('x~dx', '').replace('y~dy', '')
        self.assertNotIn('x', stripped)
        self.assertNotIn('y', stripped)


# f(x) = x: Taylor coefficients, varAt, and varOrder.
class TestLinear(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(self.x, (self.vx,), max_order=4)

    def test_at_order0(self):
        self._check(self.t.at(0), self.x)

    def test_at_order1(self):
        # f'(x)/1! = 1
        self._check(self.t.at(1), sympy.Integer(1))

    def test_at_order2_is_zero(self):
        # f''(x)/2! = 0
        self._check(self.t.at(2), sympy.Integer(0))

    def test_coeffs_length(self):
        self.assertEqual(len(self.t.coeffs), 5)

    def test_coeffs_indexed_by_order(self):
        for order in range(5):
            self._check(self.t.coeffs[(order,)], self.t.at(order))

    def test_varAt_order1_is_zero(self):
        # full_moment=moment(1)=0 and split_moments all involve moment(1)=0
        self._check(self.t.varAt(1), sympy.Integer(0))

    def test_varAt_order2(self):
        # only nn=1 contributes (d2 of x is 0 at nn=0,2)
        self._check(self.t.varAt(2), self.vx.moment(2) * self.dx**2)

    def test_varAt_order3_is_zero(self):
        # higher derivatives of x are 0; nn+pn=3 forces one >= 2
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4_is_zero(self):
        # nn+pn=4 forces one >= 2 → coeff is 0
        self._check(self.t.varAt(4), sympy.Integer(0))

    def test_varOrder_order2(self):
        self._check(self.t.varOrder(2), self.vx.moment(2) * self.dx**2)

    def test_biasAt_order1_is_zero(self):
        # biasAt(1) = δx · 1 · moment(1) = 0 (moment(1)=0)
        self._check(self.t.biasAt(1), sympy.Integer(0))

    def test_biasAt_order2_is_zero(self):
        # coeffs[(2,)] = 0 for f=x
        self._check(self.t.biasAt(2), sympy.Integer(0))

    def test_biasOrder_order2_is_zero(self):
        self._check(self.t.biasOrder(2), sympy.Integer(0))


# f(x) = x**2: Taylor coefficients, varAt, and varOrder.
class TestQuadratic(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(self.x**2, (self.vx,), max_order=4)

    def test_at_order0(self):
        self._check(self.t.at(0), self.x**2)

    def test_at_order1(self):
        # f'(x)/1! = 2x
        self._check(self.t.at(1), 2 * self.x)

    def test_at_order2(self):
        # f''(x)/2! = 1
        self._check(self.t.at(2), sympy.Integer(1))

    def test_coeffs_length(self):
        self.assertEqual(len(self.t.coeffs), 5)

    def test_coeffs_indexed_by_order(self):
        for order in range(5):
            self._check(self.t.coeffs[(order,)], self.t.at(order))

    def test_varAt_order1_is_zero(self):
        # moment(1)=0 for symmetric distribution
        self._check(self.t.varAt(1), sympy.Integer(0))

    def test_varAt_order2(self):
        # impl: j=0: x²·ζ(2)·(1-ζ(0)); j=1: 4x²·ζ(2); j=2: x²·ζ(2)·(1-ζ(0))
        m0 = self.vx.moment(0)
        m2 = self.vx.moment(2)
        x2 = self.x**2
        self._check(self.t.varAt(2),
                    self.dx**2 * (x2 * m2 * (1 - m0) + 4 * x2 * m2
                                  + x2 * m2 * (1 - m0)))

    def test_varAt_order4(self):
        # j=0,4: c_0·c_4=0 (c_4=0); j=1,3: c_1·c_3=0; j=2: 1·(ζ(4)-ζ(2)²)
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4), (m4 - m2**2) * self.dx**4)

    def test_varOrder_order1_is_zero(self):
        self._check(self.t.varOrder(1), sympy.Integer(0))

    def test_varOrder_order2(self):
        m0 = self.vx.moment(0)
        m2 = self.vx.moment(2)
        x2 = self.x**2
        self._check(self.t.varOrder(2),
                    self.dx**2 * (x2 * m2 * (1 - m0) + 4 * x2 * m2
                                  + x2 * m2 * (1 - m0)))

    def test_varOrder_order4(self):
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varOrder(4), (m4 - m2**2) * self.dx**4)

    def test_biasAt_order1_is_zero(self):
        # biasAt(1) = δx · 2x · moment(1) = 0
        self._check(self.t.biasAt(1), sympy.Integer(0))

    def test_biasAt_order2(self):
        # biasAt(2) = δx^2 · 1 · moment(2)
        self._check(self.t.biasAt(2), self.vx.moment(2) * self.dx**2)

    def test_biasAt_order3_is_zero(self):
        # coeffs[(3,)] = 0 for f=x^2
        self._check(self.t.biasAt(3), sympy.Integer(0))

    def test_biasAt_order4_is_zero(self):
        self._check(self.t.biasAt(4), sympy.Integer(0))

    def test_biasOrder_order2(self):
        self._check(self.t.biasOrder(2), self.vx.moment(2) * self.dx**2)

    # Formula (2.20)/(2.21) for f = x^c with c=2: binom(2, n) selects only n ≤ 2,
    # so the power-mean and power-precision sums collapse to short closed forms.

    def test_at_matches_formula_2_20_with_c_2(self):
        # coeffs[(n,)] = C(2, n) · x^(2-n) — the (2.20) Taylor structure
        for n in range(5):
            self._check(self.t.at(n),
                        sympy.binomial(2, n) * self.x**(2 - n))

    def test_biasAt_matches_formula_2_20_with_c_2(self):
        # biasAt(n) = δx^n · C(2, n) · x^(2-n) · ζ(n) per (2.20)
        for n in range(1, 5):
            expected = (self.dx**n * sympy.binomial(2, n)
                        * self.x**(2 - n) * self.vx.moment(n))
            self._check(self.t.biasAt(n), expected)

    def test_total_variance_matches_formula_2_21_with_c_2(self):
        # (2.21) at c=2 collapses to two terms (n=2 j=1, and n=4 j=2):
        #   δ²(x²) = 4·ζ(2)·x²·δx² + (ζ(4) - ζ(2)²)·δx⁴
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        impl_total = sum((self.t.varOrder(n) for n in range(1, 5)),
                         sympy.Integer(0))
        self._check(impl_total,
                    4 * m2 * self.x**2 * self.dx**2
                    + (m4 - m2**2) * self.dx**4)


# exp(x): n-th derivative is exp(x), so at(n) = exp(x) / n!
class TestExp(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(sympy.exp(self.x), (self.vx,), max_order=4)

    def test_order0(self):
        self._check(self.t.at(0), sympy.exp(self.x))

    def test_order1(self):
        self._check(self.t.at(1), sympy.exp(self.x))

    def test_order2(self):
        self._check(self.t.at(2), sympy.exp(self.x) / 2)

    def test_order3(self):
        self._check(self.t.at(3), sympy.exp(self.x) / 6)

    def test_varAt_order1_is_zero(self):
        self._check_varAt_order1_zero()

    def test_varAt_order2(self):
        # (2.15) at n=2, j=1: δx² · e^(2x) · ζ(2)
        self._check(self.t.varAt(2),
                    self.dx**2 * sympy.exp(self.x)**2 * self.vx.moment(2))

    def test_varAt_order3_is_zero(self):
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.15) at n=4: δx⁴ · e^(2x) · [m4/6 + (m4-m2²)/4 + m4/6]
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4),
                    self.dx**4 * sympy.exp(self.x)**2
                    * (m4 / 6 + (m4 - m2**2) / 4 + m4 / 6))

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order2(self):
        # (2.14) at n=2: δx² · (e^x/2) · ζ(2)
        self._check(self.t.biasAt(2),
                    self.dx**2 * sympy.exp(self.x) * self.vx.moment(2) / 2)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order4(self):
        # (2.14) at n=4: δx⁴ · (e^x/24) · ζ(4)
        self._check(self.t.biasAt(4),
                    self.dx**4 * sympy.exp(self.x) * self.vx.moment(4) / 24)

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)

    def test_orders_match_at_4(self):
        self._check_orders_match_at(4)


# log(x): at(0)=log(x), at(n)=(-1)^(n-1)/(n*x^n) for n>=1
class TestLog(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(sympy.log(self.x), (self.vx,), max_order=4)

    def test_order0(self):
        self._check(self.t.at(0), sympy.log(self.x))

    def test_order1(self):
        self._check(self.t.at(1), 1 / self.x)

    def test_order2(self):
        self._check(self.t.at(2), sympy.Rational(-1, 2) / self.x**2)

    def test_order3(self):
        self._check(self.t.at(3), sympy.Rational(1, 3) / self.x**3)

    def test_varAt_order1_is_zero(self):
        self._check_varAt_order1_zero()

    def test_varAt_order2(self):
        # (2.17) at n=2, j=1: P(x)²·ζ(2)/(1·1) = δx²·ζ(2)/x²
        self._check(self.t.varAt(2),
                    self.dx**2 * self.vx.moment(2) / self.x**2)

    def test_varAt_order3_is_zero(self):
        # ζ(3)=0 and ζ(j)·ζ(3-j)=0 for j∈{1,2}
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.17) at n=4: P(x)⁴·[m4/3 + (m4-m2²)/4 + m4/3]
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4),
                    self.dx**4 / self.x**4
                    * (m4 / 3 + (m4 - m2**2) / 4 + m4 / 3))

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order2(self):
        # (2.16) at n=2: P(x)² · (-1)³/2 · ζ(2) = -δx²·ζ(2)/(2x²)
        self._check(self.t.biasAt(2),
                    -self.dx**2 * self.vx.moment(2) / (2 * self.x**2))

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order4(self):
        # (2.16) at n=4: P(x)⁴ · (-1)⁵/4 · ζ(4) = -δx⁴·ζ(4)/(4x⁴)
        self._check(self.t.biasAt(4),
                    -self.dx**4 * self.vx.moment(4) / (4 * self.x**4))

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)

    def test_orders_match_at_4(self):
        self._check_orders_match_at(4)


# sin(x): derivatives cycle sin->cos->-sin->-cos->sin, divided by n!
class TestSine(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(sympy.sin(self.x), (self.vx,), max_order=4)

    def test_order0(self):
        self._check(self.t.at(0), sympy.sin(self.x))

    def test_order1(self):
        self._check(self.t.at(1), sympy.cos(self.x))

    def test_order2(self):
        self._check(self.t.at(2), -sympy.sin(self.x) / 2)

    def test_order3(self):
        self._check(self.t.at(3), -sympy.cos(self.x) / 6)

    def test_order4(self):
        self._check(self.t.at(4), sympy.sin(self.x) / 24)

    def test_varAt_order1_is_zero(self):
        self._check_varAt_order1_zero()

    def test_varAt_order2(self):
        # (2.19) at n=2, j=1: cos²(x)·ζ(2)
        self._check(self.t.varAt(2),
                    self.dx**2 * sympy.cos(self.x)**2 * self.vx.moment(2))

    def test_varAt_order3_is_zero(self):
        # ζ(3)=ζ(1)=0 so all (j∈{1,2}) terms vanish
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.19) at n=4: j=1 cos·(-cos)/6·ζ(4); j=2 sin²/4·(ζ(4)-ζ(2)²); j=3 (-cos)·cos/6·ζ(4)
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4),
                    self.dx**4 * (-sympy.cos(self.x)**2 * m4 / 3
                                  + sympy.sin(self.x)**2 * (m4 - m2**2) / 4))

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order2(self):
        # (2.18) at n=2: δx² · sin^{(2)}(x)/2! · ζ(2) = -δx²·sin(x)·ζ(2)/2
        self._check(self.t.biasAt(2),
                    -self.dx**2 * sympy.sin(self.x) * self.vx.moment(2) / 2)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order4(self):
        # (2.18) at n=4: δx⁴ · sin^{(4)}(x)/4! · ζ(4) = δx⁴·sin(x)·ζ(4)/24
        self._check(self.t.biasAt(4),
                    self.dx**4 * sympy.sin(self.x) * self.vx.moment(4) / 24)

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)

    def test_orders_match_at_4(self):
        self._check_orders_match_at(4)


# x**c (c=1.5): at(n) = c*(c-1)*...*(c-n+1)/n! * x**(c-n)
class TestPow(_Base):

    def setUp(self):
        super().setUp()
        self.c = 1.5
        self.t = analytic.StatTaylor(self.x**self.c, (self.vx,), max_order=4)

    def test_order0(self):
        self._check(self.t.at(0), self.x**self.c)

    def test_order1(self):
        self._check(self.t.at(1), self.c * self.x**(self.c - 1))

    def test_order2(self):
        self._check(self.t.at(2), self.c * (self.c - 1) / 2 * self.x**(self.c - 2))

    def test_order3(self):
        self._check(self.t.at(3), self.c * (self.c - 1) * (self.c - 2) / 6 * self.x**(self.c - 3))

    def test_varAt_order1_is_zero(self):
        self._check_varAt_order1_zero()

    def test_varAt_order2(self):
        # (2.21) at n=2, j=1: C(c,1)²·ζ(2)·x^(2c-2)·δx² = c²·x^(2c-2)·ζ(2)·δx²
        self._check(self.t.varAt(2),
                    self.dx**2 * self.x**(2 * self.c - 2)
                    * self.c**2 * self.vx.moment(2))

    def test_varAt_order3_is_zero(self):
        # ζ(3)=ζ(1)=0 so all (j∈{1,2}) terms vanish
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.21) at n=4: x^(2c-4) · sum_{j=1..3} C(c,j)·C(c,4-j)·(ζ(4)-ζ(j)·ζ(4-j))
        c = self.c
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        binom_c_2 = c * (c - 1) / 2
        binom_c_3 = c * (c - 1) * (c - 2) / 6
        self._check(self.t.varAt(4),
                    self.dx**4 * self.x**(2 * c - 4)
                    * (c * binom_c_3 * m4
                       + binom_c_2**2 * (m4 - m2**2)
                       + binom_c_3 * c * m4))

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order2(self):
        # (2.20) at n=2: δx² · C(c,2) · x^(c-2) · ζ(2) = δx²·c(c-1)/2·x^(c-2)·ζ(2)
        self._check(self.t.biasAt(2),
                    self.dx**2 * self.c * (self.c - 1) / 2
                    * self.x**(self.c - 2) * self.vx.moment(2))

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order4(self):
        # (2.20) at n=4: δx⁴ · C(c,4) · x^(c-4) · ζ(4)
        c = self.c
        self._check(self.t.biasAt(4),
                    self.dx**4 * c * (c - 1) * (c - 2) * (c - 3) / 24
                    * self.x**(c - 4) * self.vx.moment(4))

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)

    def test_orders_match_at_4(self):
        self._check_orders_match_at(4)


# f(x) = sin(x)/x: 1D rational mixture of sin(x), cos(x), and powers of x.
class TestSinXdivX(_Base):

    def setUp(self):
        super().setUp()
        self.f = sympy.sin(self.x) / self.x
        self.t = analytic.StatTaylor(self.f, (self.vx,), max_order=4)

    # at: derivatives of sin(x)/x divided by n!

    def test_order0(self):
        self._check(self.t.at(0), self.f)

    def test_order1(self):
        # f'(x) = cos(x)/x - sin(x)/x²
        self._check(self.t.at(1),
                    sympy.cos(self.x) / self.x
                    - sympy.sin(self.x) / self.x**2)

    def test_order2(self):
        # f''(x)/2! = (-sin(x)/x - 2cos(x)/x² + 2sin(x)/x³)/2
        self._check(self.t.at(2),
                    (-sympy.sin(self.x) / self.x
                     - 2 * sympy.cos(self.x) / self.x**2
                     + 2 * sympy.sin(self.x) / self.x**3) / 2)

    def test_order3(self):
        # f'''(x)/3! = (-cos(x)/x + 3sin(x)/x² + 6cos(x)/x³ - 6sin(x)/x⁴)/6
        self._check(self.t.at(3),
                    (-sympy.cos(self.x) / self.x
                     + 3 * sympy.sin(self.x) / self.x**2
                     + 6 * sympy.cos(self.x) / self.x**3
                     - 6 * sympy.sin(self.x) / self.x**4) / 6)

    def test_order4(self):
        # f''''(x)/4! = (sin(x)/x + 4cos(x)/x² - 12sin(x)/x³ - 24cos(x)/x⁴ + 24sin(x)/x⁵)/24
        self._check(self.t.at(4),
                    (sympy.sin(self.x) / self.x
                     + 4 * sympy.cos(self.x) / self.x**2
                     - 12 * sympy.sin(self.x) / self.x**3
                     - 24 * sympy.cos(self.x) / self.x**4
                     + 24 * sympy.sin(self.x) / self.x**5) / 24)

    def test_coeffs_length(self):
        self.assertEqual(len(self.t.coeffs), 5)

    # biasAt: ζ(odd) = 0 for symmetric distributions

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order2(self):
        # δx² · coeffs[(2,)] · ζ(2)
        self._check(self.t.biasAt(2),
                    self.dx**2 * self.t.at(2) * self.vx.moment(2))

    def test_biasAt_order4(self):
        # δx⁴ · coeffs[(4,)] · ζ(4)
        self._check(self.t.biasAt(4),
                    self.dx**4 * self.t.at(4) * self.vx.moment(4))

    # varAt

    def test_varAt_order1_is_zero(self):
        self._check_varAt_order1_zero()

    def test_varAt_order3_is_zero(self):
        # ζ(3)=ζ(1)=0 → all (j∈{1,2}) terms vanish
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order2(self):
        # (2.7) at n=2, j=1: c_1²·ζ(2)·δx²
        self._check(self.t.varAt(2),
                    self.dx**2 * self.t.at(1)**2 * self.vx.moment(2))

    def test_varAt_order4(self):
        # (2.7) at n=4: δx⁴·(2·c_1·c_3·ζ(4) + c_2²·(ζ(4)-ζ(2)²)) since ζ(1)=ζ(3)=0
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4),
                    self.dx**4 * (2 * self.t.at(1) * self.t.at(3) * m4
                                  + self.t.at(2)**2 * (m4 - m2**2)))

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)

    def test_orders_match_at_4(self):
        self._check_orders_match_at(4)


# f(x, y) = x + y: separable 2D function; only coeffs[(1,0)] and coeffs[(0,1)] are nonzero.
class TestXaddY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(self.x + self.y, (self.vx, self.vy), max_order=4)

    def test_at_order_0_0(self):
        self._check(self.t.at(0, 0), self.x + self.y)

    def test_at_order_1_0(self):
        # ∂f/∂x = 1
        self._check(self.t.at(1, 0), sympy.Integer(1))

    def test_at_order_0_1(self):
        # ∂f/∂y = 1
        self._check(self.t.at(0, 1), sympy.Integer(1))

    def test_at_order_1_1_is_zero(self):
        # mixed partial of a sum is 0
        self._check(self.t.at(1, 1), sympy.Integer(0))

    def test_at_order_2_0_is_zero(self):
        self._check(self.t.at(2, 0), sympy.Integer(0))

    def test_coeffs_length(self):
        # 2 vars, max_order=4: sum_{k=0..4}(k+1) = 15
        self.assertEqual(len(self.t.coeffs), 15)

    def test_varAt_order_1_1_is_zero(self):
        # full_moment=ζ_x(1)·ζ_y(1)=0; nonzero coeff pairs all have split moment with ζ(1) factor
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varAt_order_2_0(self):
        # (2.10) at p=(2,0): only nn=(1,0) contributes; c_(1,0)²·ζ_x(2)·δx²
        self._check(self.t.varAt(2, 0),
                    self.dx**2 * self.vx.moment(2))

    def test_varAt_order_0_2(self):
        # (2.10) at p=(0,2): only nn=(0,1) contributes; c_(0,1)²·ζ_y(2)·δy²
        self._check(self.t.varAt(0, 2),
                    self.dy**2 * self.vy.moment(2))

    def test_varAt_order_2_2_is_zero(self):
        # all (nn,pn) pairs use coeff that is 0
        self._check(self.t.varAt(2, 2), sympy.Integer(0))

    def test_varOrder_order2(self):
        # (2.12) δ²(x+y) = ζ_x(2)·δx² + ζ_y(2)·δy²
        self._check(self.t.varOrder(2),
                    self.vx.moment(2) * self.dx**2
                    + self.vy.moment(2) * self.dy**2)

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))

    def test_biasAt_order_1_0_is_zero(self):
        # δx · 1 · ζ_x(1) · ζ_y(0) = 0
        self._check(self.t.biasAt(1, 0), sympy.Integer(0))

    def test_biasAt_order_0_1_is_zero(self):
        self._check(self.t.biasAt(0, 1), sympy.Integer(0))

    def test_biasAt_order_1_1_is_zero(self):
        # coeffs[(1,1)] = 0
        self._check(self.t.biasAt(1, 1), sympy.Integer(0))

    def test_biasAt_order_2_0_is_zero(self):
        # coeffs[(2,0)] = 0 for f=x+y
        self._check(self.t.biasAt(2, 0), sympy.Integer(0))

    def test_biasAt_order_0_2_is_zero(self):
        # coeffs[(0,2)] = 0 for f=x+y
        self._check(self.t.biasAt(0, 2), sympy.Integer(0))

    def test_biasOrder_order2_is_zero(self):
        # all (m,n) terms with m+n=2 are zero per Formula (2.9):
        # (1,1): coeff=0; (2,0)/(0,2): coeff=0
        self._check(self.t.biasOrder(2), sympy.Integer(0))

    def test_total_bias_matches_formula_2_11(self):
        # (2.11) addition mean: x ± y is exact under ζ_k(0)=1, so total bias = 0.
        total_bias = sum((self.t.biasOrder(n) for n in range(1, 5)),
                         sympy.Integer(0))
        self._check(total_bias, sympy.Integer(0))

    def test_total_variance_matches_formula_2_12(self):
        # (2.12) δ²(x+y) = ζ_x(2)·δx² + ζ_y(2)·δy²
        total_var = sum((self.t.varOrder(n) for n in range(1, 5)),
                        sympy.Integer(0))
        self._check(total_var,
                    self.vx.moment(2) * self.dx**2
                    + self.vy.moment(2) * self.dy**2)


# f(x, y) = x * y: tests Formula (2.12) bias and (2.13) variance.
# Nonzero coeffs: (0,0)=xy, (1,0)=y, (0,1)=x, (1,1)=1; higher = 0.
class TestXmulY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(self.x * self.y, (self.vx, self.vy), max_order=4)

    def test_at_order_0_0(self):
        self._check(self.t.at(0, 0), self.x * self.y)

    def test_at_order_1_0(self):
        # ∂(xy)/∂x = y
        self._check(self.t.at(1, 0), self.y)

    def test_at_order_0_1(self):
        # ∂(xy)/∂y = x
        self._check(self.t.at(0, 1), self.x)

    def test_at_order_1_1(self):
        # ∂²(xy)/∂x∂y = 1
        self._check(self.t.at(1, 1), sympy.Integer(1))

    def test_at_order_2_0_is_zero(self):
        self._check(self.t.at(2, 0), sympy.Integer(0))

    def test_at_order_0_2_is_zero(self):
        self._check(self.t.at(0, 2), sympy.Integer(0))

    def test_coeffs_length(self):
        # 2 vars, max_order=4: sum_{k=0..4}(k+1) = 15
        self.assertEqual(len(self.t.coeffs), 15)

    # (2.13) δ²(xy) = ζ_x(2)·y²·δx² + ζ_y(2)·x²·δy² + ζ_x(2)·ζ_y(2)·δx²·δy²

    def test_varAt_order_2_0(self):
        # 1st term of (2.13): δx²·y²·ζ_x(2) (only nn=(1,0) contributes)
        self._check(self.t.varAt(2, 0),
                    self.dx**2 * self.y**2 * self.vx.moment(2))

    def test_varAt_order_0_2(self):
        # 2nd term of (2.13): δy²·x²·ζ_y(2)
        self._check(self.t.varAt(0, 2),
                    self.dy**2 * self.x**2 * self.vy.moment(2))

    def test_varAt_order_1_1_is_zero(self):
        # full=ζ_x(1)·ζ_y(1)=0; all split moments include ζ(1) factor
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varAt_order_2_1_is_zero(self):
        # only nonzero coeff pair (1,0)+(1,1) has full=ζ_x(2)·ζ_y(1)=0
        self._check(self.t.varAt(2, 1), sympy.Integer(0))

    def test_varAt_order_2_2(self):
        # 3rd term of (2.13): only nn=(1,1)=pn contributes
        self._check(self.t.varAt(2, 2),
                    self.vx.moment(2) * self.vy.moment(2)
                    * self.dx**2 * self.dy**2)

    def test_varOrder_order2(self):
        # (2,0) + (1,1)=0 + (0,2)
        self._check(self.t.varOrder(2),
                    self.dx**2 * self.y**2 * self.vx.moment(2)
                    + self.dy**2 * self.x**2 * self.vy.moment(2))

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))

    def test_varOrder_order4(self):
        # only (2,2) contributes: δx²δy²·ζ_x(2)·ζ_y(2)
        self._check(self.t.varOrder(4),
                    self.dx**2 * self.dy**2
                    * self.vx.moment(2) * self.vy.moment(2))

    def test_total_variance_matches_formula_2_13(self):
        # (2.13) δ²(xy) = ζ_x(2)·y²·δx² + ζ_y(2)·x²·δy² + ζ_x(2)·ζ_y(2)·δx²·δy²
        m2x = self.vx.moment(2)
        m2y = self.vy.moment(2)
        total = sum((self.t.varOrder(n) for n in range(1, 5)), sympy.Integer(0))
        self._check(total,
                    m2x * self.y**2 * self.dx**2
                    + m2y * self.x**2 * self.dy**2
                    + m2x * m2y * self.dx**2 * self.dy**2)

    # (2.12) overline(xy)/xy = ζ_x(0)·ζ_y(0). For Uniform ζ(0)=1, bias=0.

    def test_biasAt_order_1_0_is_zero(self):
        # δx · y · ζ_x(1) · ζ_y(0) = 0
        self._check(self.t.biasAt(1, 0), sympy.Integer(0))

    def test_biasAt_order_1_1_is_zero(self):
        # coeff=1, full_moment = ζ_x(1)·ζ_y(1) = 0
        self._check(self.t.biasAt(1, 1), sympy.Integer(0))

    def test_biasAt_order_2_0_is_zero(self):
        # coeffs[(2,0)] = 0
        self._check(self.t.biasAt(2, 0), sympy.Integer(0))

    def test_biasOrder_order2_is_zero(self):
        # (2,0)/(0,2): coeff=0; (1,1): full_moment=0
        self._check(self.t.biasOrder(2), sympy.Integer(0))


# f(x, y) = x / y. Linear in x → coeffs[(i, j)] = 0 for i ≥ 2.
# Nonzero: coeffs[(0, j)] = (-1)^j · x / y^{j+1}, coeffs[(1, j)] = (-1)^j / y^{j+1}.
class TestXdivY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(self.x / self.y, (self.vx, self.vy), max_order=4)

    # at / coeffs

    def test_at_order_0_0(self):
        self._check(self.t.at(0, 0), self.x / self.y)

    def test_at_order_1_0(self):
        # ∂(x/y)/∂x = 1/y
        self._check(self.t.at(1, 0), 1 / self.y)

    def test_at_order_0_1(self):
        # ∂(x/y)/∂y = -x/y²
        self._check(self.t.at(0, 1), -self.x / self.y**2)

    def test_at_order_1_1(self):
        # ∂²(x/y)/∂x∂y = -1/y²
        self._check(self.t.at(1, 1), -1 / self.y**2)

    def test_at_order_2_0_is_zero(self):
        # f is linear in x → all i ≥ 2 derivatives vanish
        self._check(self.t.at(2, 0), sympy.Integer(0))

    def test_at_order_0_2(self):
        self._check(self.t.at(0, 2), self.x / self.y**3)

    def test_at_order_0_3(self):
        self._check(self.t.at(0, 3), -self.x / self.y**4)

    def test_at_order_0_4(self):
        self._check(self.t.at(0, 4), self.x / self.y**5)

    def test_at_order_1_2(self):
        self._check(self.t.at(1, 2), 1 / self.y**3)

    def test_at_order_1_3(self):
        self._check(self.t.at(1, 3), -1 / self.y**4)

    def test_coeffs_length(self):
        # 2 vars, max_order=4: 15 entries
        self.assertEqual(len(self.t.coeffs), 15)

    # biasAt: nonzero iff coeffs[p] ≠ 0 AND every p_k is even (since ζ(odd)=0)
    # → only (0, 2) and (0, 4) contribute for f=x/y under symmetric distributions.

    def test_biasAt_order_1_0_is_zero(self):
        # ζ_x(1) = 0
        self._check(self.t.biasAt(1, 0), sympy.Integer(0))

    def test_biasAt_order_0_1_is_zero(self):
        self._check(self.t.biasAt(0, 1), sympy.Integer(0))

    def test_biasAt_order_1_1_is_zero(self):
        # ζ_x(1)·ζ_y(1) = 0
        self._check(self.t.biasAt(1, 1), sympy.Integer(0))

    def test_biasAt_order_2_0_is_zero(self):
        # coeffs[(2,0)] = 0
        self._check(self.t.biasAt(2, 0), sympy.Integer(0))

    def test_biasAt_order_1_2_is_zero(self):
        # ζ_x(1) = 0
        self._check(self.t.biasAt(1, 2), sympy.Integer(0))

    def test_biasAt_order_0_2(self):
        # δy² · (x/y³) · ζ_y(2)
        self._check(self.t.biasAt(0, 2),
                    self.dy**2 * self.x / self.y**3 * self.vy.moment(2))

    def test_biasAt_order_0_4(self):
        # δy⁴ · (x/y⁵) · ζ_y(4)
        self._check(self.t.biasAt(0, 4),
                    self.dy**4 * self.x / self.y**5 * self.vy.moment(4))

    def test_biasOrder_order2(self):
        # only (0,2) contributes
        self._check(self.t.biasOrder(2),
                    self.dy**2 * self.x / self.y**3 * self.vy.moment(2))

    def test_biasOrder_order3_is_zero(self):
        # (3,0): coeff=0; (0,3): ζ_y(3)=0; (1,2)/(2,1): ζ_x(1)=0 or coeff=0
        self._check(self.t.biasOrder(3), sympy.Integer(0))

    def test_biasOrder_order4(self):
        # only (0,4) contributes
        self._check(self.t.biasOrder(4),
                    self.dy**4 * self.x / self.y**5 * self.vy.moment(4))

    # varAt

    def test_varAt_order_2_0(self):
        # only nn=(1,0) pair contributes: (1/y)²·ζ_x(2)
        self._check(self.t.varAt(2, 0),
                    self.dx**2 * self.vx.moment(2) / self.y**2)

    def test_varAt_order_0_2(self):
        # only nn=(0,1) pair contributes: (-x/y²)²·ζ_y(2) = x²·ζ_y(2)/y⁴
        self._check(self.t.varAt(0, 2),
                    self.dy**2 * self.x**2 * self.vy.moment(2) / self.y**4)

    def test_varAt_order_1_1_is_zero(self):
        # full = ζ_x(1)·ζ_y(1) = 0; matching splits also vanish
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varAt_order_2_1_is_zero(self):
        # only candidate pairs have nn[0]+pn[0]=2 with both ≤1 → impossible
        # OR pair has full=ζ_y(1)=0 with vanishing split
        self._check(self.t.varAt(2, 1), sympy.Integer(0))

    def test_varAt_order_2_2(self):
        # 3 contributing pairs: (1,0)+(1,2), (1,1)+(1,1), (1,2)+(1,0)
        # each contributes c_nn·c_pn·ζ_x(2)·ζ_y(2) = ζ_x(2)·ζ_y(2)/y⁴
        self._check(self.t.varAt(2, 2),
                    3 * self.dx**2 * self.dy**2
                    * self.vx.moment(2) * self.vy.moment(2) / self.y**4)

    def test_varAt_order_0_3_is_zero(self):
        # full = ζ_y(3) = 0; all splits also vanish
        self._check(self.t.varAt(0, 3), sympy.Integer(0))

    def test_varAt_order_0_4(self):
        # contributions: (0,1)+(0,3) [×2 for symmetry] gives 2·ζ_y(4),
        # (0,2)+(0,2) gives ζ_y(4)-ζ_y(2)²; common factor x²/y⁶
        m2 = self.vy.moment(2)
        m4 = self.vy.moment(4)
        self._check(self.t.varAt(0, 4),
                    self.dy**4 * self.x**2 / self.y**6
                    * (2 * m4 + (m4 - m2**2)))

    def test_varOrder_order2(self):
        # (2,0) + (1,1)=0 + (0,2)
        self._check(self.t.varOrder(2),
                    self.dx**2 * self.vx.moment(2) / self.y**2
                    + self.dy**2 * self.x**2 * self.vy.moment(2) / self.y**4)

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))

    def test_varOrder_order4(self):
        # (4,0)=(3,1)=(1,3)=0; (2,2) and (0,4) contribute
        m2 = self.vy.moment(2)
        m4 = self.vy.moment(4)
        self._check(self.t.varOrder(4),
                    3 * self.dx**2 * self.dy**2
                    * self.vx.moment(2) * self.vy.moment(2) / self.y**4
                    + self.dy**4 * self.x**2 / self.y**6
                    * (2 * m4 + (m4 - m2**2)))


# f(x, y) = x^y. Coeffs mix x^y, x^(y-k) and powers of ln(x).
class TestXpowY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(self.x**self.y, (self.vx, self.vy), max_order=4)

    # at / coeffs

    def test_at_order_0_0(self):
        self._check(self.t.at(0, 0), self.x**self.y)

    def test_at_order_1_0(self):
        # ∂(x^y)/∂x = y · x^(y-1)
        self._check(self.t.at(1, 0), self.y * self.x**(self.y - 1))

    def test_at_order_0_1(self):
        # ∂(x^y)/∂y = x^y · log(x)
        self._check(self.t.at(0, 1), self.x**self.y * sympy.log(self.x))

    def test_at_order_2_0(self):
        # ∂²(x^y)/∂x²/2! = y(y-1)·x^(y-2)/2
        self._check(self.t.at(2, 0),
                    self.y * (self.y - 1) * self.x**(self.y - 2) / 2)

    def test_at_order_0_2(self):
        # ∂²(x^y)/∂y²/2! = x^y · log(x)²/2
        self._check(self.t.at(0, 2),
                    self.x**self.y * sympy.log(self.x)**2 / 2)

    def test_at_order_1_1(self):
        # ∂²(x^y)/∂x∂y = x^(y-1) · (1 + y·log(x))
        self._check(self.t.at(1, 1),
                    self.x**(self.y - 1) * (1 + self.y * sympy.log(self.x)))

    def test_at_order_0_4(self):
        # ∂⁴(x^y)/∂y⁴/4! = x^y · log(x)⁴/24
        self._check(self.t.at(0, 4),
                    self.x**self.y * sympy.log(self.x)**4 / 24)

    def test_coeffs_length(self):
        # 2 vars, max_order=4: 15 entries
        self.assertEqual(len(self.t.coeffs), 15)

    # biasAt: ζ(odd) = 0 zeros out any p with an odd component.

    def test_biasAt_order_1_0_is_zero(self):
        self._check(self.t.biasAt(1, 0), sympy.Integer(0))

    def test_biasAt_order_0_1_is_zero(self):
        self._check(self.t.biasAt(0, 1), sympy.Integer(0))

    def test_biasAt_order_1_1_is_zero(self):
        self._check(self.t.biasAt(1, 1), sympy.Integer(0))

    def test_biasAt_order_3_0_is_zero(self):
        self._check(self.t.biasAt(3, 0), sympy.Integer(0))

    def test_biasAt_order_2_0(self):
        # δx² · y(y-1)·x^(y-2)/2 · ζ_x(2)
        self._check(self.t.biasAt(2, 0),
                    self.dx**2 * self.y * (self.y - 1)
                    * self.x**(self.y - 2) / 2 * self.vx.moment(2))

    def test_biasAt_order_0_2(self):
        # δy² · x^y · log(x)²/2 · ζ_y(2)
        self._check(self.t.biasAt(0, 2),
                    self.dy**2 * self.x**self.y
                    * sympy.log(self.x)**2 / 2 * self.vy.moment(2))

    def test_biasOrder_order2(self):
        # (2,0) + (1,1) full=0 + (0,2)
        self._check(self.t.biasOrder(2),
                    self.dx**2 * self.y * (self.y - 1)
                    * self.x**(self.y - 2) / 2 * self.vx.moment(2)
                    + self.dy**2 * self.x**self.y
                    * sympy.log(self.x)**2 / 2 * self.vy.moment(2))

    def test_biasOrder_order3_is_zero(self):
        # every multi-index with sum=3 has at least one odd component
        self._check(self.t.biasOrder(3), sympy.Integer(0))

    # varAt

    def test_varAt_order_2_0(self):
        # only nn=(1,0) pair contributes: (y·x^(y-1))²·ζ_x(2)
        self._check(self.t.varAt(2, 0),
                    self.dx**2 * self.y**2 * self.x**(2 * self.y - 2)
                    * self.vx.moment(2))

    def test_varAt_order_0_2(self):
        # only nn=(0,1) pair contributes: (x^y·log(x))²·ζ_y(2)
        self._check(self.t.varAt(0, 2),
                    self.dy**2 * self.x**(2 * self.y)
                    * sympy.log(self.x)**2 * self.vy.moment(2))

    def test_varAt_order_1_1_is_zero(self):
        # all pairs have full = ζ_x(1)·ζ_y(1) = 0 with matching split
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varOrder_order2(self):
        # (2,0) + (1,1)=0 + (0,2)
        self._check(self.t.varOrder(2),
                    self.dx**2 * self.y**2 * self.x**(2 * self.y - 2)
                    * self.vx.moment(2)
                    + self.dy**2 * self.x**(2 * self.y)
                    * sympy.log(self.x)**2 * self.vy.moment(2))

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))


# f(x, y, z) = x + y + z: separable 3D function; nonzero coeffs only at
# (0,0,0), (1,0,0), (0,1,0), (0,0,1).
class TestXaddYaddZ(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.StatTaylor(
            self.x + self.y + self.z, (self.vx, self.vy, self.vz), max_order=4)

    def test_at_order_0_0_0(self):
        self._check(self.t.at(0, 0, 0), self.x + self.y + self.z)

    def test_at_order_1_0_0(self):
        self._check(self.t.at(1, 0, 0), sympy.Integer(1))

    def test_at_order_0_1_0(self):
        self._check(self.t.at(0, 1, 0), sympy.Integer(1))

    def test_at_order_0_0_1(self):
        self._check(self.t.at(0, 0, 1), sympy.Integer(1))

    def test_at_order_1_1_0_is_zero(self):
        # mixed partial of a sum is 0
        self._check(self.t.at(1, 1, 0), sympy.Integer(0))

    def test_at_order_2_0_0_is_zero(self):
        self._check(self.t.at(2, 0, 0), sympy.Integer(0))

    def test_coeffs_length(self):
        # 3 vars, max_order=4: sum_{k=0..4} C(k+2, 2) = 1+3+6+10+15 = 35
        self.assertEqual(len(self.t.coeffs), 35)

    def test_varAt_order_2_0_0(self):
        # 1D x-only contribution: δx² · ζ_x(2) (only nn=(1,0,0) contributes)
        self._check(self.t.varAt(2, 0, 0),
                    self.dx**2 * self.vx.moment(2))

    def test_varAt_order_0_2_0(self):
        self._check(self.t.varAt(0, 2, 0),
                    self.dy**2 * self.vy.moment(2))

    def test_varAt_order_0_0_2(self):
        self._check(self.t.varAt(0, 0, 2),
                    self.dz**2 * self.vz.moment(2))

    def test_varAt_order_1_1_0_is_zero(self):
        # full_moment includes ζ(1)·ζ(1)=0; nonzero coeff pairs all have split moment
        # involving ζ(1)=0
        self._check(self.t.varAt(1, 1, 0), sympy.Integer(0))

    def test_varAt_order_1_1_1_is_zero(self):
        self._check(self.t.varAt(1, 1, 1), sympy.Integer(0))

    def test_varOrder_order2(self):
        # 3D extension of (2.12): δ²(x+y+z) = ζ_x(2)·δx² + ζ_y(2)·δy² + ζ_z(2)·δz²
        self._check(self.t.varOrder(2),
                    self.vx.moment(2) * self.dx**2
                    + self.vy.moment(2) * self.dy**2
                    + self.vz.moment(2) * self.dz**2)

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))

    def test_biasAt_order_1_0_0_is_zero(self):
        # δx · 1 · ζ_x(1)·ζ_y(0)·ζ_z(0) = 0
        self._check(self.t.biasAt(1, 0, 0), sympy.Integer(0))

    def test_biasAt_order_2_0_0_is_zero(self):
        # coeffs[(2,0,0)] = 0
        self._check(self.t.biasAt(2, 0, 0), sympy.Integer(0))

    def test_biasAt_order_1_1_0_is_zero(self):
        # coeffs[(1,1,0)] = 0
        self._check(self.t.biasAt(1, 1, 0), sympy.Integer(0))

    def test_biasOrder_order2_is_zero(self):
        self._check(self.t.biasOrder(2), sympy.Integer(0))

    def test_total_bias_matches_formula_2_11(self):
        # 3D extension of (2.11): exact mean for x+y+z, so total bias = 0.
        total_bias = sum((self.t.biasOrder(n) for n in range(1, 5)),
                         sympy.Integer(0))
        self._check(total_bias, sympy.Integer(0))

    def test_total_variance_matches_formula_2_12(self):
        # 3D extension of (2.12): δ²(x+y+z) = ζ_x(2)·δx² + ζ_y(2)·δy² + ζ_z(2)·δz²
        total_var = sum((self.t.varOrder(n) for n in range(1, 5)),
                        sympy.Integer(0))
        self._check(total_var,
                    self.vx.moment(2) * self.dx**2
                    + self.vy.moment(2) * self.dy**2
                    + self.vz.moment(2) * self.dz**2)


class TestStatMatrix(unittest.TestCase):
    """Validates analytic.StatMatrix construction with mixed InVar/value entries
    and inheritance of index/pos/subMatrix/determ/item from the base."""

    def setUp(self):
        self.m00 = sympy.Symbol('m00')
        self.dm00 = sympy.Symbol('dm00')
        self.m11 = sympy.Symbol('m11')
        self.dm11 = sympy.Symbol('dm11')
        self.v00 = analytic.InVar(self.m00, self.dm00, analytic.EDistrType.Uniform)
        self.v11 = analytic.InVar(self.m11, self.dm11, analytic.EDistrType.Uniform)

    def test_unassigned_defaults_to_symbolic_zero(self):
        # Every (row, col) not in items resolves to sympy.S.Zero (the symbolic
        # Integer(0) singleton), and is a sympy expression rather than a Python int.
        M = analytic.StatMatrix(3, {(1, 1): self.v00})
        for r in range(3):
            for c in range(3):
                if (r, c) == (1, 1):
                    continue
                cell = M.matrix[r, c]
                self.assertIs(cell, sympy.S.Zero)
                self.assertIsInstance(cell, sympy.Expr)

    def test_mixed_entries(self):
        # 2x2 matrix: diagonal InVars, off-diagonal values 7 and 0 (default).
        M = analytic.StatMatrix(2, {(0, 0): self.v00,
                                (0, 1): 7,
                                (1, 1): self.v11})
        self.assertEqual(M.N, 2)
        self.assertEqual(len(M.in_vars), 2)
        self.assertEqual(M.matrix.tolist(),
                         [[self.m00, sympy.Integer(7)],
                          [sympy.Integer(0), self.m11]])

    def test_invalid_N(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(0, {(0, 0): self.v00})
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(-1, {(0, 0): self.v00})

    def test_no_invar_rejected(self):
        # All-value matrix: no in_vars → reject (Taylor needs at least one).
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(2, {(0, 0): 1, (1, 1): 2})

    def test_items_key_out_of_range(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(2, {(0, 0): self.v00, (5, 5): self.v11})

    def test_items_must_be_dict(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(2, [(0, 0), self.v00])

    def test_determ_with_constants(self):
        # det of [[v00, 7], [0, v11]] = v00*v11 - 0 = v00*v11
        M = analytic.StatMatrix(2, {(0, 0): self.v00,
                                (0, 1): 7,
                                (1, 1): self.v11})
        det = M.determ().function
        self.assertEqual(sympy.simplify(det - self.m00 * self.m11), 0)

    def test_item_invar_returns_taylor(self):
        M = analytic.StatMatrix(2, {(0, 0): self.v00, (1, 1): self.v11})
        T = M.item((1, 1))
        self.assertIsInstance(T, analytic.StatTaylor)
        self.assertEqual(T.function, self.m11)

    def test_item_value_position_returns_constant_taylor(self):
        # (1, 0) is the default 0 value, not an InVar; item() now returns a
        # full-in_vars StatTaylor of the entry expression (here Integer(0)).
        M = analytic.StatMatrix(2, {(0, 0): self.v00, (1, 1): self.v11})
        T = M.item((1, 0))
        self.assertIsInstance(T, analytic.StatTaylor)
        self.assertEqual(T.function, sympy.Integer(0))
        self.assertEqual(T.in_vars, M.in_vars)

    def test_subMatrix_preserves_value_entries(self):
        # 3x3 with 1 InVar and 2 nonzero value entries; drop row 0 + col 0.
        v = analytic.InVar(sympy.Symbol('a'), sympy.Symbol('da'),
                           analytic.EDistrType.Uniform)
        M = analytic.StatMatrix(3, {(1, 1): v,
                                (1, 2): 5,
                                (2, 1): 3})
        sub = M.subMatrix([(0, 0)])
        self.assertEqual(sub.N, 2)
        self.assertEqual(sub.matrix.tolist(),
                         [[sympy.Symbol('a'), sympy.Integer(5)],
                          [sympy.Integer(3), sympy.Integer(0)]])

    def test_worstmatrix_is_matrix(self):
        # Inheritance check.
        self.assertTrue(issubclass(analytic.WorstMatrix, analytic.StatMatrix))
        M = analytic.WorstMatrix(2)
        self.assertIsInstance(M, analytic.StatMatrix)

    def test_explicit_in_vars_allows_no_invar_items(self):
        # When in_vars is provided, items may consist entirely of values.
        M = analytic.StatMatrix(2, {(0, 0): 5, (1, 1): 7},
                                in_vars=(self.v00,))
        self.assertEqual(M.in_vars, (self.v00,))
        self.assertEqual(M.matrix.tolist(),
                         [[sympy.Integer(5), sympy.Integer(0)],
                          [sympy.Integer(0), sympy.Integer(7)]])

    def test_explicit_in_vars_validation(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(2, {}, in_vars=[self.v00])  # not a tuple
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(2, {}, in_vars=())  # empty
        with self.assertRaises(analytic.TaylorException):
            analytic.StatMatrix(2, {}, in_vars=('not an InVar',))


class TestWorstMatrix(unittest.TestCase):
    """Validates analytic.WorstMatrix construction, index/pos round-trips,
    subMatrix row/col deletion, and determinant equivalence with Laplace
    (cofactor) expansion along an arbitrary row or column."""

    def test_construction_2x2(self):
        M = analytic.WorstMatrix(2)
        self.assertEqual(M.N, 2)
        self.assertEqual(len(M.in_vars), 4)
        self.assertEqual(M.matrix.tolist(),
                         [[sympy.Symbol('m_0_0'), sympy.Symbol('m_0_1')],
                          [sympy.Symbol('m_1_0'), sympy.Symbol('m_1_1')]])

    def test_index_pos_roundtrip(self):
        M = analytic.WorstMatrix(4)
        for r in range(M.N):
            for c in range(M.N):
                idx = M.index(r, c)
                self.assertEqual(idx, r * M.N + c)
                self.assertEqual(M.pos(idx), (r, c))

    def test_subMatrix_drops_row_and_col(self):
        M = analytic.WorstMatrix(3)
        sub = M.subMatrix([(0, 1)])
        self.assertEqual(sub.N, 2)
        self.assertEqual(sub.matrix.tolist(),
                         [[sympy.Symbol('m_1_0'), sympy.Symbol('m_1_2')],
                          [sympy.Symbol('m_2_0'), sympy.Symbol('m_2_2')]])

    def test_subMatrix_requires_square_result(self):
        M = analytic.WorstMatrix(3)
        with self.assertRaises(analytic.TaylorException):
            M.subMatrix([(0, 0), (1, 0)])

    def test_item_returns_single_var_Taylor(self):
        M = analytic.WorstMatrix(3)
        T = M.item((1, 2))
        self.assertIsInstance(T, analytic.StatTaylor)
        self.assertEqual(T.function, sympy.Symbol('m_1_2'))
        self.assertEqual(len(T.in_vars), 1)
        self.assertIs(T.in_vars[0], M.in_vars[M.index(1, 2)])

    def test_item_default_coeffs(self):
        # at(0)=symbol, at(1)=1, at(2)=0 for a single-variable identity function
        M = analytic.WorstMatrix(2)
        T = M.item((0, 1))
        sym = sympy.Symbol('m_0_1')
        self.assertEqual(T.at(0), sym)
        self.assertEqual(T.at(1), sympy.Integer(1))
        self.assertEqual(T.at(2), sympy.Integer(0))

    def test_item_bias_and_variance(self):
        # Under normalized ζ(0)=1, bias of a single entry is 0 and variance is δm²·ζ(2).
        M = analytic.WorstMatrix(2)
        T = M.item((0, 0))
        inv = T.in_vars[0]
        self.assertEqual(sympy.simplify(T.biasOrder(0) - inv.value), 0)
        self.assertEqual(T.biasOrder(1), sympy.Integer(0))
        self.assertEqual(T.biasOrder(2), sympy.Integer(0))
        self.assertEqual(sympy.simplify(T.varOrder(2)
                                        - inv.deviation**2 * inv.moment(2)), 0)

    def test_item_invalid_position(self):
        M = analytic.WorstMatrix(2)
        with self.assertRaises(analytic.TaylorException):
            M.item((0, 1, 2))
        with self.assertRaises(analytic.TaylorException):
            M.item([0, 0])
        with self.assertRaises(analytic.TaylorException):
            M.item((2, 0))

    def test_invalid_N(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.WorstMatrix(0)
        with self.assertRaises(analytic.TaylorException):
            analytic.WorstMatrix(-1)

    def _laplace_along_row(self, M, row):
        return sum(((-1) ** (row + c) * M.matrix[row, c]
                    * M.subMatrix([(row, c)]).determ().function
                    for c in range(M.N)),
                   sympy.Integer(0))

    def _laplace_along_col(self, M, col):
        return sum(((-1) ** (r + col) * M.matrix[r, col]
                    * M.subMatrix([(r, col)]).determ().function
                    for r in range(M.N)),
                   sympy.Integer(0))

    def test_determ_matches_row_laplace_2x2(self):
        M = analytic.WorstMatrix(2)
        det = M.determ().function
        self.assertEqual(sympy.simplify(det - self._laplace_along_row(M, 0)), 0)
        self.assertEqual(sympy.simplify(det - self._laplace_along_row(M, 1)), 0)

    def test_determ_matches_col_laplace_2x2(self):
        M = analytic.WorstMatrix(2)
        det = M.determ().function
        self.assertEqual(sympy.simplify(det - self._laplace_along_col(M, 0)), 0)
        self.assertEqual(sympy.simplify(det - self._laplace_along_col(M, 1)), 0)

    def test_determ_matches_row_laplace_3x3(self):
        M = analytic.WorstMatrix(3)
        det = M.determ().function
        for r in range(M.N):
            self.assertEqual(sympy.simplify(det - self._laplace_along_row(M, r)), 0)

    def test_determ_matches_col_laplace_3x3(self):
        M = analytic.WorstMatrix(3)
        det = M.determ().function
        for c in range(M.N):
            self.assertEqual(sympy.simplify(det - self._laplace_along_col(M, c)), 0)

    def test_determ_matches_row_laplace_4x4(self):
        M = analytic.WorstMatrix(4)
        det = M.determ().function
        for r in range(M.N):
            self.assertEqual(sympy.simplify(det - self._laplace_along_row(M, r)), 0)

    def test_determ_matches_col_laplace_4x4(self):
        M = analytic.WorstMatrix(4)
        det = M.determ().function
        for c in range(M.N):
            self.assertEqual(sympy.simplify(det - self._laplace_along_col(M, c)), 0)

    def _check_adjugate_identity(self, N):
        # adj(M) · M = det(M) · I  (and also M · adj(M) = det(M) · I).
        M = analytic.WorstMatrix(N)
        adj = M.adjugate()
        self.assertIsInstance(adj, analytic.StatMatrix)
        self.assertEqual(adj.N, N)
        self.assertEqual(adj.in_vars, M.in_vars)
        det = M.determ().function
        identity = det * sympy.eye(N)
        for product in (adj.matrix * M.matrix, M.matrix * adj.matrix):
            for i in range(N):
                for j in range(N):
                    self.assertEqual(sympy.simplify(product[i, j] - identity[i, j]), 0)

    def test_adjugate_identity_1x1(self):
        self._check_adjugate_identity(1)

    def test_adjugate_identity_2x2(self):
        self._check_adjugate_identity(2)

    def test_adjugate_identity_3x3(self):
        self._check_adjugate_identity(3)

    def test_adjugate_identity_4x4(self):
        self._check_adjugate_identity(4)

    def test_adjugate_2x2_explicit_form(self):
        # For [[a, b], [c, d]], adj = [[d, -b], [-c, a]].
        M = analytic.WorstMatrix(2)
        adj = M.adjugate()
        m00 = sympy.Symbol('m_0_0')
        m01 = sympy.Symbol('m_0_1')
        m10 = sympy.Symbol('m_1_0')
        m11 = sympy.Symbol('m_1_1')
        self.assertEqual(adj.matrix.tolist(), [[m11, -m01], [-m10, m00]])

    def _check_reverse_identity(self, N):
        # M · reverse(M) = I = reverse(M) · M.
        M = analytic.WorstMatrix(N)
        inv = M.reverse()
        self.assertIsInstance(inv, analytic.StatMatrix)
        self.assertEqual(inv.N, N)
        self.assertEqual(inv.in_vars, M.in_vars)
        eye = sympy.eye(N)
        for product in (M.matrix * inv.matrix, inv.matrix * M.matrix):
            for i in range(N):
                for j in range(N):
                    self.assertEqual(sympy.simplify(product[i, j] - eye[i, j]), 0)

    def test_reverse_identity_1x1(self):
        self._check_reverse_identity(1)

    def test_reverse_identity_2x2(self):
        self._check_reverse_identity(2)

    def test_reverse_identity_3x3(self):
        self._check_reverse_identity(3)

    def test_reverse_equals_adjugate_over_det(self):
        # reverse(M) entry-wise should equal adj(M) / det(M).
        M = analytic.WorstMatrix(3)
        inv = M.reverse()
        adj = M.adjugate()
        det = M.determ().function
        for i in range(M.N):
            for j in range(M.N):
                self.assertEqual(
                    sympy.simplify(inv.matrix[i, j] - adj.matrix[i, j] / det), 0)

    def _dump_reverse_matrix(self, N: int, max_order: int = 2):
        """Dump the entire reverse(WorstMatrix(N)) via StatMatrix.dump() to
        Python/Output/dump_reverse_<N>.csv for inspection, then assert
        structural invariants of the CSV file. `max_order` is passed
        explicitly because the constructor default (16) is impractical for
        the per-entry Taylor of an N×N matrix inverse over N² in_vars."""
        M = analytic.WorstMatrix(N)
        inv = M.reverse()
        path = os.path.join(OUTDIR, 'Python', 'Output',
                            f'dump_reverse_{N}.csv')
        inv.dump(path, max_order=max_order)
        self.assertTrue(os.path.isfile(path))
        with open(path, newline='') as f:
            rows = list(csv.reader(f))
        # First N² rows are per-entry function rows (single-field, quoted by csv).
        for k, (r, c) in enumerate((rr, cc) for rr in range(N) for cc in range(N)):
            self.assertEqual(len(rows[k]), 1)
            line = rows[k][0]
            self.assertTrue(line.startswith(f'function ({r}, {c}): '))
            for inv_var in M.in_vars:
                self.assertIn(f'{inv_var.value}~{inv_var.deviation}', line)
        # Header row: order, row, col, N² InVar labels, varAt, biasAt.
        header = rows[N * N]
        self.assertEqual(header[0], 'order')
        self.assertEqual(header[1], 'row')
        self.assertEqual(header[2], 'col')
        self.assertEqual(header[-2], 'varAt')
        self.assertEqual(header[-1], 'biasAt')
        self.assertEqual(len(header), 3 + N * N + 2)
        for k, inv_var in enumerate(M.in_vars):
            self.assertEqual(header[3 + k],
                             f'{inv_var.value}~{inv_var.deviation}')
        # Every data row: order = sum of per-InVar columns; row,col in range.
        for row in rows[N * N + 1:]:
            self.assertEqual(int(row[0]),
                             sum(int(row[3 + k]) for k in range(N * N)))
            self.assertTrue(0 <= int(row[1]) < N)
            self.assertTrue(0 <= int(row[2]) < N)
        # Each (row, col) entry should produce at least one row (order-0 biasAt
        # is the entry expression, which is nonzero for every reverse entry of
        # an invertible WorstMatrix).
        positions = {(int(row[1]), int(row[2])) for row in rows[N * N + 1:]}
        self.assertEqual(positions,
                         {(r, c) for r in range(N) for c in range(N)})

    def test_dump_reverse_2x2(self):
        self._dump_reverse_matrix(2)

    def test_dump_reverse_3x3(self):
        self._dump_reverse_matrix(3)

    # The next two run a per-entry Taylor expansion of the inverse with
    # max_order=8 — for an N×N matrix that is C(N²+8, 8) multi-indices per
    # entry × N² entries × sympy.diff/varAt/biasAt of a rational function.
    # N=3 is ~24310 multi-indices × 9 entries; N=4 is ~735471 × 16. They
    # write substantial dump files into Python/Output/ for inspection.
    # Marked slow; remove the decorator (or set RUN_SLOW_DUMP_TESTS=1) to run.
    _RUN_SLOW = bool(os.environ.get('RUN_SLOW_DUMP_TESTS'))

    @unittest.skipUnless(_RUN_SLOW,
                         'slow: reverse-3x3 dump with max_order=8 (set '
                         'RUN_SLOW_DUMP_TESTS=1 to enable)')
    def test_dump_reverse_3x3_max_order_8(self):
        self._dump_reverse_matrix(3, max_order=8)

    @unittest.skipUnless(_RUN_SLOW,
                         'slow: reverse-4x4 dump with max_order=8 (set '
                         'RUN_SLOW_DUMP_TESTS=1 to enable)')
    def test_dump_reverse_4x4_max_order_8(self):
        self._dump_reverse_matrix(4, max_order=8)


if __name__ == '__main__':
    unittest.main()
