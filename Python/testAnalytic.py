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
        self.vx = analytic.InVar(self.x, self.dx, analytic.EDistrType.Uniform)
        self.vy = analytic.InVar(self.y, self.dy, analytic.EDistrType.Uniform)

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

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order2(self):
        # δx² · exp(x)/2 · ζ(2, κ)
        self._check(self.t.biasAt(2),
                    self.dx**2 * sympy.exp(self.x) * self.vx.moment(2) / 2)

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)


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

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order2(self):
        # δx² · (-1/(2x²)) · ζ(2, κ)
        self._check(self.t.biasAt(2),
                    -self.dx**2 * self.vx.moment(2) / (2 * self.x**2))

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)


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

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order2(self):
        # δx² · (-sin(x)/2) · ζ(2, κ)
        self._check(self.t.biasAt(2),
                    -self.dx**2 * sympy.sin(self.x) * self.vx.moment(2) / 2)

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)


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

    def test_biasAt_order1_is_zero(self):
        self._check_biasAt_odd_zero(1)

    def test_biasAt_order3_is_zero(self):
        self._check_biasAt_odd_zero(3)

    def test_biasAt_order2(self):
        # δx² · c(c-1)/2 · x^(c-2) · ζ(2, κ)
        self._check(self.t.biasAt(2),
                    self.dx**2 * self.c * (self.c - 1) / 2
                    * self.x**(self.c - 2) * self.vx.moment(2))

    def test_orders_match_at_2(self):
        self._check_orders_match_at(2)


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

    def test_varAt_order_2_2_is_zero(self):
        # all (nn,pn) pairs use coeff that is 0
        self._check(self.t.varAt(2, 2), sympy.Integer(0))

    def test_varOrder_order2_is_zero(self):
        # only p=(1,1) contributes at order 2
        self._check(self.t.varOrder(2), sympy.Integer(0))

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
        # coeffs[(2,0)] = 0
        self._check(self.t.biasAt(2, 0), sympy.Integer(0))

    def test_biasOrder_order2_is_zero(self):
        self._check(self.t.biasOrder(2), sympy.Integer(0))

