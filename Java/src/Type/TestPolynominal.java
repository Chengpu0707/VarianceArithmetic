package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.junit.Test;

public class TestPolynominal {

    @Test
    public void test_poly_0() {
        VarDbl res;
        try {
            res = new VarDbl().polynominal(new double[] {1});
            assertEquals(res.value(), 1, 0);
            assertEquals(res.uncertainty(), 0, 0);

            res = new VarDbl(1).polynominal(new double[] {1});
            assertEquals(res.value(), 1, 1);
            assertEquals(res.uncertainty(), 0, 0);

            res = new VarDbl(1).polynominal(new double[] {-2});
            assertEquals(res.value(), -2, 0);
            assertEquals(res.uncertainty(), 0, 0);

            res = new VarDbl(2).polynominal(new double[] {1.0});
            assertEquals(res.value(), 1, 0);
            assertEquals(res.uncertainty(), 0, 0);
        } catch (InitException | DivergentException | NotReliableException e) {
            fail(e.getMessage());
        }
    }
    
    @Test
    public void test_poly_1() {
        VarDbl res;
        try {
            res = new VarDbl(0, 1/8).polynominal(new double[]{0, 1});
            assertEquals(res.value(), 0, 0);
            assertEquals(res.uncertainty(), 1/8, 1e-6);

            res = new VarDbl(0, 1/8).polynominal(new double[]{1,1});
            assertEquals(res.value(), 1, 0);
            assertEquals(res.uncertainty(), 1/8, 1e-6);

            res = new VarDbl(-2, 1/8).polynominal(new double[]{1,1});
            assertEquals(res.value(), -1, 0);
            assertEquals(res.uncertainty(), 1/8, 1e-6);
        } catch (InitException | DivergentException | NotReliableException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void test_poly_2() {
        VarDbl res, res2;
        try {
            for (double value: new double[]{0, 1, -1, 2, -2, 0.25, -0.25}) {
                res = new VarDbl(value, 0.5).polynominal(new double[]{0,0,1});
                res2 = new VarDbl(value*value + 0.5*0.5, 4 *value*value *0.5*0.5 + 2 *0.5*0.5*0.5*0.5, true);
                assertEquals(res.value(), res2.value(), 1e-5);
                assertEquals(res.uncertainty(), res2.uncertainty(), 5e-5);

                res = new VarDbl(value, 0.5).polynominal(new double[]{1,2,1});
                res2 = new VarDbl((value + 1)*(value + 1) + 0.5*0.5, 4 *(value + 1)*(value + 1) *0.5*0.5 + 2 *0.5*0.5*0.5*0.5, true);
                assertEquals(res.value(), res2.value(), 1e-5);
                assertEquals(res.uncertainty(), res2.uncertainty(), 5e-5);

                res = new VarDbl(value, 0.5).polynominal(new double[]{1,-2,1});
                res2 = new VarDbl((value - 1)*(value - 1) + 0.5*0.5, 4 *(value - 1)*(value - 1) *0.5*0.5 + 2 *0.5*0.5*0.5*0.5, true);
                assertEquals(res.value(), res2.value(), 1e-5);
                assertEquals(res.uncertainty(), res2.uncertainty(), 5e-5);
            }
        } catch (InitException | DivergentException | NotReliableException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void test_poly_3() {
        VarDbl res, res2;
        try {
            for (double value: new double[]{0, 1, -1, 2, -2, 0.25, -0.25}) {
                res = new VarDbl(value, 0.5).polynominal(new double[]{0,0,0,1});
                res2 = new VarDbl(value, 0.5).power(3);
                assertEquals(res.value(), res2.value(), 0);
                assertEquals(res.uncertainty(), res2.uncertainty(), 0);

                res = new VarDbl(value, 0.5).polynominal(new double[]{1,3,3,1});
                res2 = new VarDbl(1 + value, 0.5).power(3);
                assertEquals(res.value(), res2.value(), 0);
                assertEquals(res.uncertainty(), res2.uncertainty(), 0);

                res = new VarDbl(value, 0.5).polynominal(new double[]{1,-3,3,-1});
                res2 = new VarDbl(1 - value, 0.5).power(3);
                assertEquals(res.value(), res2.value(), 0);
                assertEquals(res.uncertainty(), res2.uncertainty(), 0);
            }
        } catch (InitException | DivergentException | NotReliableException | IOException | NotMonotonicException | NotStableException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testExpansion() {
        try {
            VarDbl res = new VarDbl(2, 0.2).polynominal(new double[]{0,0,1}, "./Java/Output/Poly_1_0.5_2.txt");
            try (BufferedReader br = new BufferedReader(new FileReader("./Java/Output/Poly_1_0.5_2.txt"))) {
                String line, last = null;
                int cnt = 0;
                while ((line = br.readLine()) != null) {
                    ++cnt;
                    last = line;
                }
                assertEquals(12, cnt);
                final String[] sLast = last.split("\t");
                assertEquals(res.value(), Double.valueOf(sLast[0]), 1e-5);          
                assertEquals(res.variance(), Double.valueOf(sLast[1]) + Double.valueOf(sLast[2]), 5e-5);
            } catch (IOException e) {
                fail(e.getMessage());
            }
        } catch (InitException | DivergentException | NotReliableException | IOException e) {
            fail(e.getMessage());
        }
    }

}
