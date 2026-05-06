"""Unit tests for analytic.py (statistical Taylor expansion).

Validation tests:
  TestInVar              — InVar construction, defaults, immutability, type checks.
  TestInVarMoment        — InVar.moment() for Gaussian (symbolic ζ) and Uniform (numeric).
  TestTaylor             — Taylor construction, properties, immutability, type checks.
  TestTaylorMethod       — Taylor.at() argument validation and small structural cases.
  TestTaylorVarAt        — varAt(*orders) argument validation and 1D/2D structural cases.
  TestTaylorVarOrder     — varOrder(n) argument validation and structural sums.
  TestTaylorBiasAt       — biasAt(*orders) argument validation (incl. all-zero rejection).
  TestTaylorBiasOrder    — biasOrder(n) argument validation.

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

import math
import sympy
import unittest

import analytic


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
        import momentum as _mmt
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertAlmostEqual(v.moment(0), _mmt.Uniform(bounding=math.sqrt(3))[0])

    def test_uniform_default_kappa_order2(self):
        import momentum as _mmt
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertAlmostEqual(v.moment(2), _mmt.Uniform(bounding=math.sqrt(3))[2])

    def test_uniform_odd_order_is_zero(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.assertEqual(v.moment(1), 0)

    def test_custom_kappa_gaussian(self):
        v = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian, kappa=3.0)
        self.assertEqual(v.moment(2), analytic.zeta(2, v.kappa))


class TestTaylor(unittest.TestCase):

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.y = sympy.Symbol('y')
        self.dy = sympy.Symbol('dy')
        self.vx = analytic.InVar(self.x, self.dx, analytic.EDistrType.Gaussian)
        self.vy = analytic.InVar(self.y, self.dy, analytic.EDistrType.Uniform)

    def test_single_var(self):
        t = analytic.Taylor(sympy.sin(self.x), (self.vx,), max_order=2)
        self.assertEqual(t.function, sympy.sin(self.x))
        self.assertEqual(t.in_vars, (self.vx,))

    def test_two_vars(self):
        t = analytic.Taylor(self.x + self.y, (self.vx, self.vy), max_order=2)
        self.assertEqual(len(t.in_vars), 2)

    def test_default_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=200)
        self.assertEqual(t.max_order, 200)

    def test_custom_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=5)
        self.assertEqual(t.max_order, 5)

    def test_readonly_function(self):
        t = analytic.Taylor(sympy.sin(self.x), (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.function = sympy.cos(self.x)

    def test_readonly_in_vars(self):
        t = analytic.Taylor(sympy.sin(self.x), (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.in_vars = (self.vy,)

    def test_readonly_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=5)
        with self.assertRaises(AttributeError):
            t.max_order = 10

    def test_readonly_coeffs(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.coeffs = (self.x,)

    def test_no_extra_attributes(self):
        t = analytic.Taylor(sympy.sin(self.x), (self.vx,), max_order=2)
        with self.assertRaises(AttributeError):
            t.extra = 42

    def test_invalid_function_type(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor('sin(x)', (self.vx,))

    def test_invalid_in_vars_not_tuple(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(sympy.sin(self.x), [self.vx])

    def test_invalid_in_vars_empty(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(sympy.sin(self.x), ())

    def test_invalid_in_vars_element(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(sympy.sin(self.x), (self.x,))

    def test_invalid_max_order_float(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(self.x, (self.vx,), max_order=2.0)

    def test_invalid_max_order_bool(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(self.x, (self.vx,), max_order=True)

    def test_invalid_max_order_negative(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(self.x, (self.vx,), max_order=-1)


class TestTaylorMethod(unittest.TestCase):

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
        t = analytic.Taylor(self.x*self.y, (self.vx, self.vy), max_order=2)
        self._check(t.at(1, 0), self.y)
        self._check(t.at(0, 1), self.x)

    def test_two_vars_order2(self):
        # f(x,y)=xy: pure partials ∂²f/∂x²=0, ∂²f/∂y²=0 (mixed term not stored)
        t = analytic.Taylor(self.x*self.y, (self.vx, self.vy), max_order=2)
        self._check(t.at(2, 0), sympy.Integer(0))
        self._check(t.at(0, 2), sympy.Integer(0))

    def test_invalid_function_contains_deviation(self):
        with self.assertRaises(analytic.TaylorException):
            analytic.Taylor(self.x + self.dx, (self.vx,), max_order=2)

    def test_taylor_exceeds_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(3)

    def test_invalid_order_float(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(1.0)

    def test_invalid_order_bool(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(True)

    def test_invalid_order_negative(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=2)
        with self.assertRaises(analytic.TaylorException):
            t.at(-1)


class _Base(unittest.TestCase):
    # Uniform default kappa=√3 gives moment(0)=1 and moment(2)=1, clean symbolic results.

    def setUp(self):
        self.x = sympy.Symbol('x')
        self.dx = sympy.Symbol('dx')
        self.y = sympy.Symbol('y')
        self.dy = sympy.Symbol('dy')
        self.z = sympy.Symbol('z')
        self.dz = sympy.Symbol('dz')
        self.vx = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.vy = analytic.InVar(self.y, self.dy, analytic.EDistrType.Uniform)
        self.vz = analytic.InVar(self.z, self.dz, analytic.EDistrType.Uniform)

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


class TestTaylorVarAt(_Base):

    def test_invalid_wrong_num_orders(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(1, 1)

    def test_invalid_order_float(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(2.0)

    def test_invalid_order_bool(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(True)

    def test_invalid_order_zero(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(0)

    def test_invalid_order_exceeds_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varAt(5)

    def test_2d_product_order2_2(self):
        # f(x,y)=x*y: p=(2,2), nn=(1,1) contributes (mixed partials beyond order 1 are 0)
        t = analytic.Taylor(self.x * self.y, (self.vx, self.vy), max_order=4)
        self._check(t.varAt(2, 2),
                    self.vx.moment(2) * self.vy.moment(2) * self.dx**2 * self.dy**2)


class TestTaylorVarOrder(_Base):

    def test_invalid_n_float(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(2.0)

    def test_invalid_n_bool(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(True)

    def test_invalid_n_zero(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(0)

    def test_invalid_n_exceeds_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.varOrder(5)

    def test_2d_product_order4(self):
        # f(x,y)=x*y: only varAt(2,2) contributes among (1,3),(2,2),(3,1)
        t = analytic.Taylor(self.x * self.y, (self.vx, self.vy), max_order=4)
        self._check(t.varOrder(4),
                    self.vx.moment(2) * self.vy.moment(2) * self.dx**2 * self.dy**2)


class TestTaylorBiasAt(_Base):

    def test_invalid_wrong_num_orders(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(1, 1)

    def test_invalid_order_float(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(2.0)

    def test_invalid_order_bool(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(True)

    def test_invalid_order_negative(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(-1)

    def test_invalid_all_zero(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(0)

    def test_invalid_2d_all_zero(self):
        t = analytic.Taylor(self.x + self.y, (self.vx, self.vy), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(0, 0)

    def test_invalid_order_exceeds_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasAt(5)

    def test_2d_partial_zero_allowed(self):
        # biasAt(1,0)=δx · ∂f/∂x · ζ_x(1) · ζ_y(0); for f=x and symmetric x, ζ_x(1)=0
        t = analytic.Taylor(self.x, (self.vx, self.vy), max_order=4)
        self._check(t.biasAt(1, 0), sympy.Integer(0))


class TestTaylorBiasOrder(_Base):

    def test_invalid_n_float(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(2.0)

    def test_invalid_n_bool(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(True)

    def test_invalid_n_zero(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(0)

    def test_invalid_n_exceeds_max_order(self):
        t = analytic.Taylor(self.x, (self.vx,), max_order=4)
        with self.assertRaises(analytic.TaylorException):
            t.biasOrder(5)


# f(x) = x: Taylor coefficients, varAt, and varOrder.
class TestLinear(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.Taylor(self.x, (self.vx,), max_order=4)

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
        self.t = analytic.Taylor(self.x**2, (self.vx,), max_order=4)

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
        # nn=0 and nn=2 cancel (moment(0)=1 for Uniform)
        self._check(self.t.varAt(2), 4 * self.x**2 * self.vx.moment(2) * self.dx**2)

    def test_varAt_order4(self):
        # only nn=2 contributes (higher derivatives of x^2 are 0)
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4), (m4 - m2**2) * self.dx**4)

    def test_varOrder_order1_is_zero(self):
        self._check(self.t.varOrder(1), sympy.Integer(0))

    def test_varOrder_order2(self):
        self._check(self.t.varOrder(2), 4 * self.x**2 * self.vx.moment(2) * self.dx**2)

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


# exp(x): n-th derivative is exp(x), so at(n) = exp(x) / n!
class TestExp(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.Taylor(sympy.exp(self.x), (self.vx,), max_order=4)

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
        # (2.15) at n=2: δx² · (e^x)² · (ζ(2) - ζ(1)²)/(1!·1!) = δx² · (e^x)² · ζ(2)
        # Equals our impl when ζ(0)=1 (Uniform default).
        self._check(self.t.varAt(2),
                    self.dx**2 * sympy.exp(self.x)**2 * self.vx.moment(2))

    def test_varAt_order3_is_zero(self):
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.15) at n=4: δx⁴ · (e^x)² · sum_{j=1..3} (ζ(4) - ζ(j)·ζ(4-j))/(j!(4-j)!)
        # j=1,3: each (ζ(4) - 0)/6 since ζ(1)=ζ(3)=0; j=2: (ζ(4) - ζ(2)²)/4
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
        self.t = analytic.Taylor(sympy.log(self.x), (self.vx,), max_order=4)

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
        # (2.17) at n=2: P(x)² · (ζ(2) - ζ(1)²)/(1·1) = δx²·ζ(2)/x²
        self._check(self.t.varAt(2),
                    self.dx**2 * self.vx.moment(2) / self.x**2)

    def test_varAt_order3_is_zero(self):
        # ζ(3)=0 and ζ(j)·ζ(3-j)=0 for j∈{1,2}
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.17) at n=4: P(x)⁴ · sum_{j=1..3} (ζ(4) - ζ(j)·ζ(4-j))/(j·(4-j))
        # j=1,3: ζ(4)/3 each (ζ(1)=ζ(3)=0); j=2: (ζ(4) - ζ(2)²)/4
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
        self.t = analytic.Taylor(sympy.sin(self.x), (self.vx,), max_order=4)

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
        # (2.19) at n=2, j=1: cos(x)² · (ζ(2) - ζ(1)²)/(1!·1!) = cos²(x)·ζ(2)
        self._check(self.t.varAt(2),
                    self.dx**2 * sympy.cos(self.x)**2 * self.vx.moment(2))

    def test_varAt_order3_is_zero(self):
        # ζ(3)=ζ(1)=0 so all (j∈{1,2}) terms vanish
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.19) at n=4: sum_{j=1..3} sin^{(j)}(x)·sin^{(4-j)}(x) · (ζ(4) - ζ(j)·ζ(4-j))/(j!·(4-j)!)
        # j=1: cos·(-cos)/6 · ζ(4); j=2: sin²/4 · (ζ(4) - ζ(2)²); j=3: (-cos)·cos/6 · ζ(4)
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4),
                    self.dx**4 * (
                        -sympy.cos(self.x)**2 * m4 / 3
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
        self.t = analytic.Taylor(self.x**self.c, (self.vx,), max_order=4)

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
        # (2.21) at n=2, j=1: C(c,1)² · (ζ(2) - ζ(1)²) = c²·ζ(2); times x^(2c-2)
        self._check(self.t.varAt(2),
                    self.dx**2 * self.x**(2 * self.c - 2)
                    * self.c**2 * self.vx.moment(2))

    def test_varAt_order3_is_zero(self):
        # ζ(3)=ζ(1)=0 so all (j∈{1,2}) terms vanish
        self._check(self.t.varAt(3), sympy.Integer(0))

    def test_varAt_order4(self):
        # (2.21) at n=4: sum_{j=1..3} C(c,j)·C(c,4-j) · (ζ(4) - ζ(j)·ζ(4-j))
        # j=1: c·c(c-1)(c-2)/6 · ζ(4); j=2: [c(c-1)/2]² · (ζ(4) - ζ(2)²);
        # j=3: c(c-1)(c-2)/6·c · ζ(4); times x^(2c-4)
        c = self.c
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        self._check(self.t.varAt(4),
                    self.dx**4 * self.x**(2 * c - 4)
                    * (c**2 * (c - 1) * (c - 2) * m4 / 6
                       + c**2 * (c - 1)**2 * (m4 - m2**2) / 4
                       + c**2 * (c - 1) * (c - 2) * m4 / 6))

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
        self.t = analytic.Taylor(self.f, (self.vx,), max_order=4)

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
        # only j=1 contributes: coeffs[(1,)]² · (ζ(2) - ζ(1)²) = coeffs[(1,)]² · ζ(2)
        # Equals our impl when ζ(0)=1 (Uniform default).
        self._check(self.t.varAt(2),
                    self.dx**2 * self.t.at(1)**2 * self.vx.moment(2))

    def test_varAt_order4(self):
        # j=1,3: coeffs[(1,)]·coeffs[(3,)]·ζ(4); j=2: coeffs[(2,)]²·(ζ(4)-ζ(2)²)
        m2 = self.vx.moment(2)
        m4 = self.vx.moment(4)
        c1 = self.t.at(1)
        c2 = self.t.at(2)
        c3 = self.t.at(3)
        self._check(self.t.varAt(4),
                    self.dx**4 * (c1 * c3 * m4
                                  + c2**2 * (m4 - m2**2)
                                  + c3 * c1 * m4))

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)

    def test_orders_match_at_4(self):
        self._check_orders_match_at(4)


# f(x, y) = x + y: separable 2D function; only coeffs[(1,0)] and coeffs[(0,1)] are nonzero.
class TestXaddY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.Taylor(self.x + self.y, (self.vx, self.vy), max_order=4)

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
        # full_moment=ζ_x(1)·ζ_y(1)=0; nonzero coeff pairs all have split moment with ζ(1) factor
        self._check(self.t.varAt(2, 0), self.vx.moment(2) * self.dx**2)

    def test_varAt_order_0_2(self):
        # full_moment=ζ_x(1)·ζ_y(1)=0; nonzero coeff pairs all have split moment with ζ(1) factor
        self._check(self.t.varAt(0, 2), self.vy.moment(2) * self.dy**2)

    def test_varAt_order_2_2_is_zero(self):
        # all (nn,pn) pairs use coeff that is 0
        self._check(self.t.varAt(2, 2), sympy.Integer(0))

    def test_varOrder_order2(self):
        # only p=(1,1) contributes at order 2
        self._check(self.t.varOrder(2),  self.vx.moment(2) * self.dx**2 + self.vy.moment(2) * self.dy**2)

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


# f(x, y) = x * y: tests Formula (2.12) bias and (2.13) variance.
# Nonzero coeffs: (0,0)=xy, (1,0)=y, (0,1)=x, (1,1)=1; higher = 0.
class TestXmulY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.Taylor(self.x * self.y, (self.vx, self.vy), max_order=4)

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
        # 1st term of (2.13): nn=(1,0)=pn contributes y²·ζ_x(2)
        self._check(self.t.varAt(2, 0),
                    self.vx.moment(2) * self.y**2 * self.dx**2)

    def test_varAt_order_0_2(self):
        # 2nd term of (2.13)
        self._check(self.t.varAt(0, 2),
                    self.vy.moment(2) * self.x**2 * self.dy**2)

    def test_varAt_order_1_1_is_zero(self):
        # full=ζ_x(1)·ζ_y(1)=0; all split moments include ζ(1) factor
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varAt_order_2_1_is_zero(self):
        # only nonzero coeff pair (1,0)+(1,1) has full=ζ_x(2)·ζ_y(1)=0
        # and split also vanishes
        self._check(self.t.varAt(2, 1), sympy.Integer(0))

    def test_varAt_order_2_2(self):
        # 3rd term of (2.13): only nn=(1,1)=pn contributes
        # coeff² · (ζ_x(2)·ζ_y(2) - 0) = ζ_x(2)·ζ_y(2)
        self._check(self.t.varAt(2, 2),
                    self.vx.moment(2) * self.vy.moment(2)
                    * self.dx**2 * self.dy**2)

    def test_varOrder_order2(self):
        # 1st + 2nd terms of (2.13): varAt(2,0) + varAt(1,1) + varAt(0,2)
        self._check(self.t.varOrder(2),
                    self.vx.moment(2) * self.y**2 * self.dx**2
                    + self.vy.moment(2) * self.x**2 * self.dy**2)

    def test_varOrder_order3_is_zero(self):
        # (3,0),(0,3): coeffs all 0; (2,1),(1,2): see test_varAt_order_2_1
        self._check(self.t.varOrder(3), sympy.Integer(0))

    def test_varOrder_order4(self):
        # 3rd term of (2.13): only varAt(2,2) contributes among (m,n) with m+n=4
        self._check(self.t.varOrder(4),
                    self.vx.moment(2) * self.vy.moment(2)
                    * self.dx**2 * self.dy**2)

    def test_total_variance_matches_2_13(self):
        # full Formula (2.13): δ²(xy) = sum of varOrder(n) over all n
        total = sum((self.t.varOrder(n) for n in range(1, 5)), sympy.Integer(0))
        m2x = self.vx.moment(2)
        m2y = self.vy.moment(2)
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
        self.t = analytic.Taylor(self.x / self.y, (self.vx, self.vy), max_order=4)

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
        # δy² · (x/y³) · ζ_x(0) · ζ_y(2)
        self._check(self.t.biasAt(0, 2),
                    self.dy**2 * self.x * self.vy.moment(2) / self.y**3)

    def test_biasAt_order_0_4(self):
        # δy⁴ · (x/y⁵) · ζ_x(0) · ζ_y(4)
        self._check(self.t.biasAt(0, 4),
                    self.dy**4 * self.x * self.vy.moment(4) / self.y**5)

    def test_biasOrder_order2(self):
        # only (0,2) contributes
        self._check(self.t.biasOrder(2),
                    self.dy**2 * self.x * self.vy.moment(2) / self.y**3)

    def test_biasOrder_order3_is_zero(self):
        # (3,0): coeff=0; (0,3): ζ_y(3)=0; (1,2)/(2,1): ζ_x(1)=0 or coeff=0
        self._check(self.t.biasOrder(3), sympy.Integer(0))

    def test_biasOrder_order4(self):
        # only (0,4) contributes among (m,n) with m+n=4
        self._check(self.t.biasOrder(4),
                    self.dy**4 * self.x * self.vy.moment(4) / self.y**5)

    # varAt

    def test_varAt_order_2_0(self):
        # only nn=(1,0)=pn contributes: (1/y)² · ζ_x(2) = ζ_x(2)/y²
        self._check(self.t.varAt(2, 0),
                    self.vx.moment(2) * self.dx**2 / self.y**2)

    def test_varAt_order_0_2(self):
        # only nn=(0,1)=pn contributes: (-x/y²)² · ζ_y(2) = x²·ζ_y(2)/y⁴
        self._check(self.t.varAt(0, 2),
                    self.vy.moment(2) * self.x**2 * self.dy**2 / self.y**4)

    def test_varAt_order_1_1_is_zero(self):
        # full = ζ_x(1)·ζ_y(1) = 0; matching splits also vanish
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varAt_order_2_1_is_zero(self):
        # only candidate pairs have nn[0]+pn[0]=2 with both ≤1 → impossible
        # OR pair has full=ζ_y(1)=0 with vanishing split
        self._check(self.t.varAt(2, 1), sympy.Integer(0))

    def test_varAt_order_2_2(self):
        # only nn[0]=pn[0]=1 pairs contribute. Three j-splits each give m2x·m2y/y⁴.
        m2x = self.vx.moment(2)
        m2y = self.vy.moment(2)
        self._check(self.t.varAt(2, 2),
                    3 * m2x * m2y * self.dx**2 * self.dy**2 / self.y**4)

    def test_varAt_order_0_3_is_zero(self):
        # full = ζ_y(3) = 0; all splits also vanish
        self._check(self.t.varAt(0, 3), sympy.Integer(0))

    def test_varAt_order_0_4(self):
        # nn=(0,1)+(0,3): x²·m4y/y⁶; nn=(0,2)+(0,2): x²·(m4y-m2y²)/y⁶;
        # nn=(0,3)+(0,1): x²·m4y/y⁶ (term order matches impl's accumulation)
        m2y = self.vy.moment(2)
        m4y = self.vy.moment(4)
        common = self.dy**4 * self.x**2 / self.y**6
        self._check(self.t.varAt(0, 4),
                    common * m4y + common * (m4y - m2y**2) + common * m4y)

    def test_varOrder_order2(self):
        # varAt(2,0) + varAt(1,1) + varAt(0,2)
        self._check(self.t.varOrder(2),
                    self.vx.moment(2) * self.dx**2 / self.y**2
                    + self.vy.moment(2) * self.x**2 * self.dy**2 / self.y**4)

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))

    def test_varOrder_order4(self):
        # varAt(2,2) + varAt(0,4); (4,0)/(3,1)/(1,3) all 0
        m2x = self.vx.moment(2)
        m2y = self.vy.moment(2)
        m4y = self.vy.moment(4)
        common = self.dy**4 * self.x**2 / self.y**6
        self._check(self.t.varOrder(4),
                    3 * m2x * m2y * self.dx**2 * self.dy**2 / self.y**4
                    + common * m4y + common * (m4y - m2y**2) + common * m4y)


# f(x, y) = x^y. Coeffs mix x^y, x^(y-k) and powers of ln(x).
class TestXpowY(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.Taylor(self.x**self.y, (self.vx, self.vy), max_order=4)

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
        # only (2,0) and (0,2) contribute (others have ζ(1)=0)
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
        # nn=(1,0)=pn: coeffs² = (y·x^(y-1))² and (full - split) = ζ_x(2)
        self._check(self.t.varAt(2, 0),
                    self.dx**2 * self.y**2
                    * self.x**(2 * self.y - 2) * self.vx.moment(2))

    def test_varAt_order_0_2(self):
        # nn=(0,1)=pn: coeffs² = (x^y·log(x))² and (full - split) = ζ_y(2)
        self._check(self.t.varAt(0, 2),
                    self.dy**2 * self.x**(2 * self.y)
                    * sympy.log(self.x)**2 * self.vy.moment(2))

    def test_varAt_order_1_1_is_zero(self):
        # all pairs have full = ζ_x(1)·ζ_y(1) = 0 with matching split
        self._check(self.t.varAt(1, 1), sympy.Integer(0))

    def test_varOrder_order2(self):
        # only (2,0) and (0,2) contribute
        self._check(self.t.varOrder(2),
                    self.dx**2 * self.y**2
                    * self.x**(2 * self.y - 2) * self.vx.moment(2)
                    + self.dy**2 * self.x**(2 * self.y)
                    * sympy.log(self.x)**2 * self.vy.moment(2))

    def test_varOrder_order3_is_zero(self):
        self._check(self.t.varOrder(3), sympy.Integer(0))


# f(x, y, z) = x + y + z: separable 3D function; nonzero coeffs only at
# (0,0,0), (1,0,0), (0,1,0), (0,0,1).
class TestXaddYaddZ(_Base):

    def setUp(self):
        super().setUp()
        self.t = analytic.Taylor(
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
        # 1D x-only contribution: δx² · ζ_x(2)
        self._check(self.t.varAt(2, 0, 0), self.vx.moment(2) * self.dx**2)

    def test_varAt_order_0_2_0(self):
        self._check(self.t.varAt(0, 2, 0), self.vy.moment(2) * self.dy**2)

    def test_varAt_order_0_0_2(self):
        self._check(self.t.varAt(0, 0, 2), self.vz.moment(2) * self.dz**2)

    def test_varAt_order_1_1_0_is_zero(self):
        # full_moment includes ζ(1)·ζ(1)=0; nonzero coeff pairs all have split moment
        # involving ζ(1)=0
        self._check(self.t.varAt(1, 1, 0), sympy.Integer(0))

    def test_varAt_order_1_1_1_is_zero(self):
        self._check(self.t.varAt(1, 1, 1), sympy.Integer(0))

    def test_varOrder_order2(self):
        # only (2,0,0), (0,2,0), (0,0,2) contribute; cross-pair multi-indices give 0
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


if __name__ == '__main__':
    unittest.main()
