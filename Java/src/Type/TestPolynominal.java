package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.IOException;

import org.junit.Test;


public class TestPolynominal {

    @Test
    public void test_poly_0() 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, InitException {
        VarDbl res, res2;
        for (double value: new double[]{0, 1, -1, 2, -2, 0.25, -0.25}) {
            res2 = new VarDbl(value);
            res = res2.polynominal(new double[] {-2});
            assertEquals(res.value(), -2, 0);
            assertEquals(res.uncertainty(), 0, 0);
            res = new VarDbl(-2).polynominal(new VarDbl[] {res2});
            assertEquals(res.value(), res2.value(), 0);
            assertEquals(res.uncertainty(), res2.uncertainty(), 0);
        }
    }
    
    @Test
    public void test_poly_1() 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, InitException {
        VarDbl res, res2;
        for (double value: new double[]{0, 1, -1, 2, -2, 0.25, -0.25}) {
            res = new VarDbl(value).polynominal(new double[] {0, -2});
            res2 = new VarDbl(value).multiply(-2);
            assertEquals(res.value(), res2.value(), 0);
            assertEquals(res.uncertainty(), res2.uncertainty(), 0);
            assertEquals(res.value(), -2 *value, 0);
            assertEquals(res.uncertainty(), 0, 0);

            res = new VarDbl(-2).polynominal(new VarDbl[] {new VarDbl(), new VarDbl(value)});
            res2 = new VarDbl(-2).multiply(new VarDbl(value));
            assertEquals(res.value(), res2.value(), 0);
            assertEquals(res.uncertainty(), res2.uncertainty(), 0);
            assertEquals(res.value(), -2 *value, 0);
            assertEquals(res.uncertainty(), 0, 0);

            res = new VarDbl(value, 0.5).polynominal(new double[] {0, -2});
            res2 = new VarDbl(value, 0.5).multiply(-2);
            assertEquals(res.value(), -2*value, 0);
            assertEquals(res.uncertainty(), 1, 1e-5);

            res = new VarDbl(-2).polynominal(new VarDbl[] {new VarDbl(), new VarDbl(value, 0.5)});
            res2 = new VarDbl(-2).multiply(new VarDbl(value, 0.5));
            assertEquals(res.value(), res2.value(), 0);
            assertEquals(res.uncertainty(), res2.uncertainty(), 0);
            assertEquals(res.value(), -2 *value, 0);
            assertEquals(res.uncertainty(), 1, 1e-5);

        }
    }

    @Test
    public void test_poly_2() 
            throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                    InitException, IOException {
        VarDbl res, res2;
        for (double value: new double[]{0, 1, -1, 2, -2, 0.25, -0.25}) {
            res = new VarDbl(value, 0.5).polynominal(new double[]{0,0,1});
            res2 = new VarDbl(value, 0.5).pow(2);
            assertEquals(res.value(), res2.value(), 1e-5);
            assertEquals(res.uncertainty(), res2.uncertainty(), 5e-5);

            res = new VarDbl(value, 0.5).polynominal(new double[]{1,2,1});
            res2 = new VarDbl(value + 1, 0.5).pow(2);
            assertEquals(res.value(), res2.value(), 1e-5);
            assertEquals(res.uncertainty(), res2.uncertainty(), 5e-5);

            res = new VarDbl(value, 0.5).polynominal(new double[]{1,-2,1});
            res2 = new VarDbl(value - 1, 0.5).pow(2);
            assertEquals(res.value(), res2.value(), 1e-5);
            assertEquals(res.uncertainty(), res2.uncertainty(), 5e-5);
        }
    }

    @Test
    public void test_poly_3() {
        VarDbl res, res2;
        try {
            for (double value: new double[]{0, 1, -1, 2, -2, 0.25, -0.25}) {
                res = new VarDbl(value, 0.5).polynominal(new double[]{0,0,0,1});
                res2 = new VarDbl(value, 0.5).pow(3);
                assertEquals(res.value(), res2.value(), 0);
                assertEquals(res.uncertainty(), res2.uncertainty(), 0);

                res = new VarDbl(value, 0.5).polynominal(new double[]{1,3,3,1});
                res2 = new VarDbl(1 + value, 0.5).pow(3);
                assertEquals(res.value(), res2.value(), 0);
                assertEquals(res.uncertainty(), res2.uncertainty(), 0);

                res = new VarDbl(value, 0.5).polynominal(new double[]{1,-3,3,-1});
                res2 = new VarDbl(1 - value, 0.5).pow(3);
                assertEquals(res.value(), res2.value(), 0);
                assertEquals(res.uncertainty(), res2.uncertainty(), 0);
            }
        } catch (Throwable e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testExpansion() 
        throws NotFiniteException, NotReliableException, NotMonotonicException, NotStableException, NotPositiveException, 
                InitException, IOException {
        final String dumpFile = "./Java/Output/Poly_1_0.5_2.txt";
        final double[] s1dPoly = new double[]{1,-2,1};
        VarDbl res1 = new VarDbl(2, 0.5).polynominal(s1dPoly, dumpFile);
        VarDbl res2 = new VarDbl(1, 0.5).pow(2);
        assertEquals(res1.value(), res2.value(), 0);
        assertEquals(res1.uncertainty(), res2.uncertainty(), 0);

        VarDbl.DumpResult res = VarDbl.readDumpFile(dumpFile, res1, null);
        assertEquals(res.s1dPoly.sDbl.length, s1dPoly.length);
        for (int i = 0; i < s1dPoly.length; ++i) {
            assertEquals(res.s1dPoly.sDbl[i], s1dPoly[i], 0);
        }
    }

}
