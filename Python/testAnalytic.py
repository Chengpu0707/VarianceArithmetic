import sympy
import unittest

import analytic

class TestSeries (unittest.TestCase):
    '''
    The series is Maclaurin series
    '''
    def testSin(self):
        x = sympy.symbols('x', is_real=True)
        sin = sympy.sin(x)
        self.assertTupleEqual(sin.args, (x,))
        e = sin.series(x, 0, 10).removeO()
        self.assertDictEqual(e.as_coefficients_dict(), {
            x: sympy.Integer(1), x**3: sympy.Rational(-1, 6), 
            x**5: sympy.Rational(1, 120), x**7: sympy.Rational(-1, 5040), 
            x**9: sympy.Rational(1, 362880)
        })


    def testPow(self):
        x, c = sympy.symbols('x c')
        pow = sympy.Pow(x, c)
        self.assertTupleEqual(pow.args, (x, c))
        e = pow.series(x, 0, 10).removeO()
        self.assertDictEqual(e.as_coefficients_dict(), {
            x**c: 1,
        })

        pow = sympy.Pow(x + 1, c)
        self.assertTupleEqual(pow.args, (x + 1, c))
        e = pow.series(x, 0, 5).removeO()
        self.assertDictEqual(e.as_coefficients_dict(), {
            1: 1,
            x *c: 1,
            x**2 *c*(c-1): sympy.Rational(1, 2),
            x**3 *c*(c-1)*(c-2): sympy.Rational(1, 6),
            x**4 *c*(c-1)*(c-2)*(c-3): sympy.Rational(1, 24),
        })


class TestCompare (unittest.TestCase):

    def testNameOrder(self):
        '''
        x < y by default
        '''
        x, y = sympy.symbols('x y', real=True)
        self.assertEqual(x.compare(y), -1)
        self.assertEqual(y.compare(x), 1)
        with self.assertRaises(TypeError):
            bool(x < y)

    def testAssumption(self):
        '''
        Assumption can only be True
        '''
        x = sympy.Symbol('x', positive=True)
        self.assertTrue(0 < x)
        y = sympy.Symbol('y', negative=True)
        self.assertTrue(y < 0)
        self.assertTrue(y < x)
        self.assertTrue(y <= x)
        z = sympy.Symbol('z', positive=False)
        with self.assertRaises(TypeError):
            bool(z <= 0)
        z = sympy.Symbol('z', nonpositive=True)
        self.assertTrue(z <= 0)

    def testInequality(self):
        x = sympy.symbols('x', real=True)
        self.assertEqual(sympy.Pow(x, 2).compare(x), 1)
        self.assertEqual(x.compare(sympy.Pow(x, 2)), -1)


class TestTaylor1d (unittest.TestCase):
    MAX_ORDER = 9

    def _testBiasVar(self, sDiff, sBias, sVar):
        self.assertEqual(sympy.simplify(sBias[(2,)] - sDiff[(2,)] *analytic.momentum(2)), 0)
        self.assertEqual(sympy.simplify(sBias[(4,)] - sDiff[(4,)] *analytic.momentum(4)), 0)
        self.assertEqual(sympy.simplify(sBias[(6,)] - sDiff[(6,)] *analytic.momentum(6)), 0)
        self.assertEqual(sympy.simplify(sBias[(8,)] - sDiff[(8,)] *analytic.momentum(8)), 0)

        eq = sDiff[(1,)]**2 *analytic.momentum(2)
        self.assertEqual(sympy.simplify(sVar[(2,)] - eq), 0)
        eq = (sDiff[(2,)]**2 + sDiff[(1,)]*sDiff[(3,)]*2) *analytic.momentum(4) - sBias[(2,)]**2
        self.assertEqual(sympy.simplify(sVar[(4,)] - eq), 0)
        eq = (sDiff[(3,)]**2 + sDiff[(2,)]*sDiff[(4,)]*2 + sDiff[(1,)]*sDiff[(5,)]*2) *analytic.momentum(6) \
             - sBias[(2,)]*sBias[(4,)]*2
        self.assertEqual(sympy.simplify(sVar[(6,)] - eq), 0)
        eq = (sDiff[(4,)]**2 + sDiff[(3,)]*sDiff[(5,)]*2 + sDiff[(2,)]*sDiff[(6,)]*2 + sDiff[(1,)]*sDiff[(7,)]*2) *analytic.momentum(8) \
             - sBias[(2,)]*sBias[(6,)]*2 - sBias[(4,)]**2
        self.assertEqual(sympy.simplify(sVar[(8,)] - eq), 0)


    def testExp(self):
        x = sympy.symbols('x', is_real=True)
        exp = sympy.exp(x)
        sDiff, sBias, sVar = analytic.taylor_series(exp, (x,), maxOrder=TestTaylor1d.MAX_ORDER)

        self.assertEqual(sDiff[(0,)], exp)
        self.assertEqual(sDiff[(1,)], exp)
        self.assertEqual(sDiff[(2,)], exp/2)
        self.assertEqual(sDiff[(3,)], exp/6)
        self.assertEqual(sDiff[(4,)], exp/24)
        self.assertEqual(sDiff[(5,)], exp/120)
        self.assertEqual(sDiff[(6,)], exp/720)
        self.assertEqual(sDiff[(7,)], exp/5040)
        self.assertEqual(sDiff[(8,)], exp/40320)

        self._testBiasVar(sDiff, sBias, sVar)

    def testSin(self):
        x = sympy.symbols('x', is_real=True)
        sin = sympy.sin(x)
        cos = sympy.cos(x)
        sDiff, sBias, sVar = analytic.taylor_series(sin, (x,), maxOrder=TestTaylor1d.MAX_ORDER)

        self.assertEqual(sDiff[(0,)], sin)
        self.assertEqual(sDiff[(1,)], +cos)
        self.assertEqual(sDiff[(2,)], -sin/2)
        self.assertEqual(sDiff[(3,)], -cos/6)
        self.assertEqual(sDiff[(4,)], +sin/24)
        self.assertEqual(sDiff[(5,)], +cos/120)
        self.assertEqual(sDiff[(6,)], -sin/720)
        self.assertEqual(sDiff[(7,)], -cos/5040)
        self.assertEqual(sDiff[(8,)], +sin/40320)

        self._testBiasVar(sDiff, sBias, sVar)

    def testLog(self):
        x = sympy.symbols('x', is_positive=True)
        log = sympy.log(x)
        sDiff, sBias, sVar = analytic.taylor_series(log, (x,), maxOrder=TestTaylor1d.MAX_ORDER)

        self.assertEqual(sDiff[(0,)], log)
        self.assertEqual(sDiff[(1,)], +1/x)
        self.assertEqual(sDiff[(2,)], -sympy.Rational(1, 2) /x**2)
        self.assertEqual(sDiff[(3,)], +sympy.Rational(1, 3) /x**3)
        self.assertEqual(sDiff[(4,)], -sympy.Rational(1, 4) /x**4)
        self.assertEqual(sDiff[(5,)], +sympy.Rational(1, 5) /x**5)
        self.assertEqual(sDiff[(6,)], -sympy.Rational(1, 6) /x**6)
        self.assertEqual(sDiff[(7,)], +sympy.Rational(1, 7) /x**7)
        self.assertEqual(sDiff[(8,)], -sympy.Rational(1, 8) /x**8)

        self._testBiasVar(sDiff, sBias, sVar)

    def testPow(self):
        x = sympy.symbols('x', is_positive=True)
        c = sympy.symbols('c', is_real=True)
        pow = sympy.Pow(x, c)
        sDiff, sBias, sVar = analytic.taylor_series(pow, (x,), maxOrder=TestTaylor1d.MAX_ORDER)

        self.assertEqual(sDiff[(0,)], pow)
        self.assertEqual(sympy.simplify(sDiff[(1,)] - pow /x**1 * c), 0)
        self.assertEqual(sympy.simplify(sDiff[(2,)] - pow /x**2 * c*(c-1) /2), 0)
        self.assertEqual(sympy.simplify(sDiff[(3,)] - pow /x**3 * c*(c-1)*(c-2) /6), 0)
        self.assertEqual(sympy.simplify(sDiff[(4,)] - pow /x**4 * c*(c-1)*(c-2)*(c-3) /24), 0)
        self.assertEqual(sympy.simplify(sDiff[(5,)] - pow /x**5 * c*(c-1)*(c-2)*(c-3)*(c-4) /120), 0)
        self.assertEqual(sympy.simplify(sDiff[(6,)] - pow /x**6 * c*(c-1)*(c-2)*(c-3)*(c-4)*(c-5) /720), 0)
        self.assertEqual(sympy.simplify(sDiff[(7,)] - pow /x**7 * c*(c-1)*(c-2)*(c-3)*(c-4)*(c-5)*(c-6) /5040), 0)
        self.assertEqual(sympy.simplify(sDiff[(8,)] - pow /x**8 * c*(c-1)*(c-2)*(c-3)*(c-4)*(c-5)*(c-6)*(c-7) /40320), 0)

        self._testBiasVar(sDiff, sBias, sVar)


class TestTaylor2d (unittest.TestCase):

    def testAdd(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x + y
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, {(0,0): func, (1,0): 1, (0,1): 1})
        self.assertDictEqual(sBias, {})
        self.assertDictEqual(sVar, {(2,0): analytic.momentum(2), (0,2): analytic.momentum(2)})

    def testMinus(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x - y
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, {(0,0): func, (1,0): 1, (0,1): -1})
        self.assertDictEqual(sBias, {})
        self.assertDictEqual(sVar, {(2,0): analytic.momentum(2), (0,2): analytic.momentum(2)})

    def testMul(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x * y
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, {(0,0): func, (1,0): y, (0,1): x, (1,1): sympy.Integer(1)})
        self.assertDictEqual(sBias, {})
        self.assertSetEqual(set(sVar.keys()), set([(2,0), (0,2), (2,2)]))
        self.assertEqual(sympy.simplify(sVar[(2,0)] - y**2 *analytic.momentum(2)), 0)
        self.assertEqual(sympy.simplify(sVar[(0,2)] - x**2 *analytic.momentum(2)), 0)
        self.assertEqual(sympy.simplify(sVar[(2,2)] - analytic.momentum(2)**2), 0)

    def testX2(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x**2
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, {(0,0): func, (1,0): x*2, (2,0): sympy.Integer(1)})
        self.assertDictEqual(sBias, {(2,0): analytic.momentum(2)})
        self.assertSetEqual(set(sVar.keys()), set([(2,0), (4,0)]))
        eq = x**2 * 4 *analytic.momentum(2)
        self.assertEqual(sympy.simplify(sVar[(2,0)] - eq), 0)
        eq = analytic.momentum(4) - analytic.momentum(2)**2
        self.assertEqual(sympy.simplify(sVar[(4,0)] - eq), 0)

    def testY2(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = y**2
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, {(0,0): func, (0,1): y*2, (0,2): sympy.Integer(1)})
        self.assertDictEqual(sBias, {(0,2): analytic.momentum(2)})
        self.assertSetEqual(set(sVar.keys()), set([(0,2), (0,4)]))
        eq = y**2 * 4 *analytic.momentum(2)
        self.assertEqual(sympy.simplify(sVar[(0,2)] - eq), 0)
        eq = analytic.momentum(4) - analytic.momentum(2)**2
        self.assertEqual(sympy.simplify(sVar[(0,4)] - eq), 0)

    def testX2Y2(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x**2 + y**2
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, { (0,0): func, 
            (1,0): x*2, (2,0): sympy.Integer(1),
            (0,1): y*2, (0,2): sympy.Integer(1)
        })
        self.assertDictEqual(sBias, {(2,0): analytic.momentum(2), (0,2): analytic.momentum(2)})
        self.assertSetEqual(set(sVar.keys()), set([(2,0), (0,2), (4,0), (0,4), (2,2)]))
        self.assertEqual(sympy.simplify(sVar[(2,0)] - x**2 * 4 *analytic.momentum(2)), 0)
        self.assertEqual(sympy.simplify(sVar[(0,2)] - y**2 * 4 *analytic.momentum(2)), 0)
        self.assertEqual(sympy.simplify(sVar[(4,0)] -(analytic.momentum(4) - analytic.momentum(2)**2)), 0)
        self.assertEqual(sympy.simplify(sVar[(0,4)] -(analytic.momentum(4) - analytic.momentum(2)**2)), 0)
        self.assertEqual(sympy.simplify(sVar[(2,2)]), 0)

    def testSum2(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x**2 + y**2 + 2*x*y
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, { (0,0): func, 
            (1,0): (x + y)*2, (2,0): sympy.Integer(1),
            (0,1): (x + y)*2, (0,2): sympy.Integer(1),
            (1,1): sympy.Integer(2)
        }) 
        self.assertDictEqual(sBias, {(2,0): analytic.momentum(2), (0,2): analytic.momentum(2)})
        self.assertSetEqual(set(sVar.keys()), set([(2,0), (0,2), (4,0), (0,4), (2,2)]))
        eq = (x + y)**2 * 4 *analytic.momentum(2)
        self.assertEqual(sympy.simplify(sVar[(2,0)] - eq), 0)
        self.assertEqual(sympy.simplify(sVar[(0,2)] - eq), 0)
        eq = analytic.momentum(4) - analytic.momentum(2)**2
        self.assertEqual(sympy.simplify(sVar[(4,0)] - eq), 0)
        self.assertEqual(sympy.simplify(sVar[(0,4)] - eq), 0)
        eq = 4 *analytic.momentum(2)**2
        self.assertEqual(sympy.simplify(sVar[(2,2)] - eq), 0)

    def testMinus2(self):
        x, y = sympy.symbols('x y', is_real=True)
        func = x**2 + y**2 - 2*x*y
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=7)
        self.assertDictEqual(sDiff, { (0,0): func, 
            (1,0): (x - y)*2, (2,0): sympy.Integer(1),
            (0,1): (y - x)*2, (0,2): sympy.Integer(1),
            (1,1): sympy.Integer(-2)
        }) 
        self.assertDictEqual(sBias, {(2,0): analytic.momentum(2), (0,2): analytic.momentum(2)})
        self.assertSetEqual(set(sVar.keys()), set([(2,0), (0,2), (4,0), (0,4), (2,2)]))
        eq = (x - y)**2 * 4 *analytic.momentum(2)
        self.assertEqual(sympy.simplify(sVar[(2,0)] - eq), 0)
        self.assertEqual(sympy.simplify(sVar[(0,2)] - eq), 0)
        eq = analytic.momentum(4) - analytic.momentum(2)**2
        self.assertEqual(sympy.simplify(sVar[(4,0)] - eq), 0)
        self.assertEqual(sympy.simplify(sVar[(0,4)] - eq), 0)
        eq = 4 *analytic.momentum(2)**2
        self.assertEqual(sympy.simplify(sVar[(2,2)] - eq), 0)

    def testPower(self):
        x, y = sympy.symbols('x y', is_positive=True)
        func = x**y
        sDiff, sBias, sVar = analytic.taylor_series(func, (x, y), maxOrder=5)
        self.assertDictEqual(sDiff, {
            (0, 0): x**y, 
            (1, 0): x**y *y/x, 
            (0, 1): x**y *sympy.log(x), 
            (2, 0): x**y *y**2/(2*x**2) - x**y *y/(2*x**2), 
            (1, 1): x**y *y*sympy.log(x)/x + x**y /x,
            (0, 2): x**y *sympy.log(x)**2 /2, 
            (3, 0): x**y*y**3/(6*x**3) - x**y*y**2/(2*x**3) + x**y*y/(3*x**3), 
            (2, 1): x**y*y**2*sympy.log(x)/(2*x**2) - x**y*y*sympy.log(x)/(2*x**2) + x**y*y/x**2 - x**y/(2*x**2), 
            (1, 2): x**y*y*sympy.log(x)**2/(2*x) + x**y*sympy.log(x)/x, 
            (0, 3): x**y*sympy.log(x)**3/6, 
            (4, 0): x**y*y**4/(24*x**4) - x**y*y**3/(4*x**4) + 11*x**y*y**2/(24*x**4) - x**y*y/(4*x**4), 
            (3, 1): x**y*y**3*sympy.log(x)/(6*x**3) - x**y*y**2*sympy.log(x)/(2*x**3) + x**y*y**2/(2*x**3) 
                    + x**y*y*sympy.log(x)/(3*x**3) - x**y*y/x**3 + x**y/(3*x**3), 
            (2, 2): x**y*y**2*sympy.log(x)**2/(4*x**2) - x**y*y*sympy.log(x)**2/(4*x**2) 
                    + x**y*y*sympy.log(x)/x**2 - x**y*sympy.log(x)/(2*x**2) + x**y/(2*x**2), 
            (1, 3): x**y*y*sympy.log(x)**3/(6*x) + x**y*sympy.log(x)**2/(2*x), 
            (0, 4): x**y*sympy.log(x)**4/24, 
        })
        self.assertDictEqual(sBias, {
            (2, 0): (x**y *y**2 /(2*x**2) - x**y *y/(2*x**2)) *analytic.momentum(2), 
            (0, 2): x**y *sympy.log(x)**2 *analytic.momentum(2)/2, 
            (4, 0): (x**y*y**4/(24*x**4) - x**y*y**3/(4*x**4) + 11*x**y*y**2/(24*x**4) 
                     - x**y*y/(4*x**4))*analytic.momentum(4), 
            (2, 2): (x**y*y**2*sympy.log(x)**2/(4*x**2) - x**y*y*sympy.log(x)**2/(4*x**2) 
                     + x**y*y*sympy.log(x)/x**2 - x**y*sympy.log(x)/(2*x**2) + x**y/(2*x**2))*analytic.momentum(4), 
            (0, 4): x**y*sympy.log(x)**4*analytic.momentum(4)/24
        })
        self.assertDictEqual(sVar, {
            (2, 0): x**(2*y) *y**2 *analytic.momentum(2)/x**2, 
            (0, 2): x**(2*y) *sympy.log(x)**2 *analytic.momentum(2), 
            (4, 0): -(x**y*y**2/(2*x**2) - x**y*y/(2*x**2))**2*analytic.momentum(2)**2 
                    + (x**y*y**2/(2*x**2) - x**y*y/(2*x**2))**2*analytic.momentum(4) 
                    + 2*x**y*y*(x**y*y**3/(6*x**3) - x**y*y**2/(2*x**3) + x**y*y/(3*x**3))*analytic.momentum(4)/x, 
            (2, 2): 2*x**y*(x**y*y**2*sympy.log(x)/(2*x**2) - x**y*y*sympy.log(x)/(2*x**2) 
                    + x**y*y/x**2 - x**y/(2*x**2))*sympy.log(x)*analytic.momentum(2)**2 
                    + (x**y*y*sympy.log(x)/x + x**y/x)**2*analytic.momentum(2)**2 
                    + 2*x**y*y*(x**y*y*sympy.log(x)**2/(2*x) + x**y*sympy.log(x)/x)*analytic.momentum(2)**2/x,
            (0, 4): -x**(2*y) *sympy.log(x)**4 *analytic.momentum(2)**2 /4 
                    + 7* x**(2*y)* sympy.log(x)**4 *analytic.momentum(4)/12,
        })
        return



class TestMatrix (unittest.TestCase):
    def testMatrixSymbol(self):
        '''
        Do not how to order {m[0,0], m[0,1], m[1,0], m[1,1]}
        Do not know how to expand det as m[0,0]*m[1,1] - m[0,1]*m[1,0]
        Expand det with [m[0,0], m[0,1], m[1,0], m[1,1] will cause infinitive looping
        '''
        m = sympy.MatrixSymbol('m', 2, 2)
        self.assertSetEqual(m.free_symbols, {m})
        self.assertEqual(m.as_explicit(), sympy.Matrix([[m[0,0], m[0,1]], [m[1,0], m[1,1]]]))
        self.assertSetEqual(m.as_explicit().free_symbols, {m})
        self.assertSetEqual(set(m), {m[0,0], m[0,1], m[1,0], m[1,1]})
#        self.assertListEqual(list(set(m)), [m[0,0], m[0,1], m[1,0], m[1,1]])
        det = m.det()
        self.assertNotEqual(sympy.simplify(det - (m[0,0]*m[1,1] - m[0,1]*m[1,0])), 0)
        diff = sympy.diff(det, m[0,0])
        self.assertEqual(sympy.simplify(diff - sympy.Derivative(sympy.Determinant(m), m[0,0])), 0)



    def test2x2(self):
        w, x, y, z = sympy.symbols('w x y z')
        mat = sympy.Matrix([[w, x], [y, z]])
        mat.free_symbols

        det = mat.det()
        sDiff, sBias, sVar = analytic.taylor_series(det, (w, x, y, z))
        self.assertEqual(sympy.simplify(det - (w*z - x*y)), 0)
        self.assertDictEqual(sDiff, {
            (0,0,0,0): w*z - x*y,
            (1,0,0,0): z, (0,1,0,0): -y, (0,0,1,0): -x, (0,0,0,1): w, 
            (1,0,0,1): 1, (0,1,1,0): -1
        })
        self.assertDictEqual(sBias, {})
        self.assertDictEqual(sVar, { 
            (2,0,0,0): z**2 *analytic.momentum(2), (0,2,0,0): y**2 *analytic.momentum(2), 
            (0,0,2,0): x**2 *analytic.momentum(2), (0,0,0,2): w**2 *analytic.momentum(2),
            (2,0,0,2): analytic.momentum(2)**2, (0,2,2,0): analytic.momentum(2)**2
        })

        adj = mat.adjugate()
        self.assertListEqual(adj.tolist(), [[z, -x], [-y, w]])
        sDiff, sBias, sVar = analytic.taylor_series(adj, (w, x, y, z))
        self.assertDictEqual(dict(sDiff[0,0]), {(0, 0, 0, 0):  z, (0, 0, 0, 1):  1})
        self.assertDictEqual(dict(sDiff[0,1]), {(0, 0, 0, 0): -x, (0, 1, 0, 0): -1})
        self.assertDictEqual(dict(sDiff[1,0]), {(0, 0, 0, 0): -y, (0, 0, 1, 0): -1})
        self.assertDictEqual(dict(sDiff[1,1]), {(0, 0, 0, 0):  w, (1, 0, 0, 0):  1})
        self.assertDictEqual(dict(sBias[0,0]), {})
        self.assertDictEqual(dict(sBias[0,1]), {})
        self.assertDictEqual(dict(sBias[1,0]), {})
        self.assertDictEqual(dict(sBias[1,1]), {})
        self.assertDictEqual(dict(sVar[0,0]), {(0, 0, 0, 2): analytic.momentum(2)})
        self.assertDictEqual(dict(sVar[0,1]), {(0, 2, 0, 0): analytic.momentum(2)})
        self.assertDictEqual(dict(sVar[1,0]), {(0, 0, 2, 0): analytic.momentum(2)})
        self.assertDictEqual(dict(sVar[1,1]), {(2, 0, 0, 0): analytic.momentum(2)})
 




if __name__ == '__main__':
    unittest.main()