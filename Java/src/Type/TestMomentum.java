
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;


public class TestMomentum {

    /*
     * https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
     */

    @Test
    public void testCdf() {
        assertEquals(0, Momentum.cdf(-9), 0);
        assertEquals(0.0013498980316302145, Momentum.cdf(-3), 2E-16);
        assertEquals(0.02275013194817932, Momentum.cdf(-2), 2E-16);
        assertEquals(0.15865525393145702, Momentum.cdf(-1), 2E-16);
        assertEquals(0.5, Momentum.cdf(0), 0);
        assertEquals(0.84134474606854298, Momentum.cdf(1), 2E-16);      // 0.68268949213708596
        assertEquals(0.97724986805182068, Momentum.cdf(2), 2E-16);      // 0.95449973610364136
        assertEquals(0.9986501019683697855, Momentum.cdf(3), 2E-16);    // 0.997300203936739571
        assertEquals(1, Momentum.cdf(9), 0);

        assertEquals(0.8140331597529114, Momentum.cdf(0.002/0.00224), 2E-16);
    }

    @Test
    public void testInverseCdf() {
        assertEquals(-3, Momentum.inverseCDF(0.0013498980316302145), 1E-8);
        assertEquals(-2, Momentum.inverseCDF(0.02275013194817932), 1E-8);
        assertEquals(-1, Momentum.inverseCDF(0.15865525393145702), 1E-8);
        assertEquals(0, Momentum.inverseCDF(0.5), 1E-8);
        assertEquals(+1, Momentum.inverseCDF(0.84134474606854298), 1E-8);
        assertEquals(+2, Momentum.inverseCDF(0.97724986805182068), 1E-8);
        assertEquals(+3, Momentum.inverseCDF(0.9986501019683697855), 1E-8);
    }

    @Test
    public void TestFactorial() {
        assertEquals(1, Momentum.factorial(0), 0);
        assertEquals(1, Momentum.factorial(1), 0);
        assertEquals(2, Momentum.factorial(2), 0);
        assertEquals(6, Momentum.factorial(3), 0);
        assertEquals(24, Momentum.factorial(4), 0);
        assertEquals(120, Momentum.factorial(5), 0);
        assertEquals(9.33262154439441E155, Momentum.factorial(99), 1E155);
        assertEquals(7.257415615307994E306, Momentum.factorial(170), 1E306);
        assertEquals(Double.POSITIVE_INFINITY, Momentum.factorial(170)*171, 0);
    }

    @Test
    public void TestDoubleFactorial() {
        assertEquals(1, Momentum.doubleFactorial(0), 0);
        assertEquals(1, Momentum.doubleFactorial(1), 0);
        assertEquals(2, Momentum.doubleFactorial(2), 0);
        assertEquals(3, Momentum.doubleFactorial(3), 0);
        assertEquals(8, Momentum.doubleFactorial(4), 0);
        assertEquals(15, Momentum.doubleFactorial(5), 0);
        assertEquals(2.7253921397507295E78, Momentum.doubleFactorial(99), 1E78);
        assertEquals(1.0898143681335199E154, Momentum.doubleFactorial(170), 1E154);
        assertEquals(1.863582569508319E156, Momentum.doubleFactorial(170)*171, 1E156);
    }

    @Test
    public void TestResidual24() {
        double x = 3;
        int n = 24;
        final boolean odd = ((n % 2) == 1);
        double term = odd?  1 : x;
        term *= Momentum.pdf(x);
        double sum = term;
        assertEquals(0.0044318484119380075, Momentum.pdf(x), 2E-16);
        assertEquals(0.013295545235814023, term, 2E-16);
         /*
        for (int j = odd? 2 : 3; j < n; j += 2) {
            term *= x * x / j;
            sum += term;
        }
         */
        int j = 3;
        term *= x * x / j;
        sum += term;
        assertEquals(0.03988663570744207, term, 2E-16);
        j = 5;
        term *= x * x / j;
        sum += term;
        assertEquals(0.07179594427339572, term, 2E-16);
        j = 7;
        term *= x * x / j;
        sum += term;
        assertEquals(0.09230907120865164, term, 2E-16);
        j = 9;
        term *= x * x / j;
        sum += term;
        assertEquals(0.09230907120865164, term, 2E-16);
        assertEquals(0.3095962676339551, sum, 2E-16);
        double sum2 = 0;
        for (j = 11; j < n; j += 2) {
            term *= x * x / j;
            sum2 += term;
        }
        assertEquals(0.1883519332126673, sum2, 2E-16);
        assertEquals(0.4979482008466224, sum + sum2, 2E-16);

        assertEquals(0.4979482008466224, Momentum.residual(n, x), 2E-16);
    }

    @Test
    public void TestResidual25() {
        double x = 3;
        int n = 25;
        final boolean odd = ((n % 2) == 1);
        double term = odd?  1 : x;
        term *= Momentum.pdf(x);
        double sum = term;
        assertEquals(0.0044318484119380075, Momentum.pdf(x), 2E-16);
        assertEquals(0.0044318484119380075, term, 2E-16);
         /*
        for (int j = odd? 2 : 3; j < n; j += 2) {
            term *= x * x / j;
            sum += term;
        }
         */
        int j = 2;
        term *= x * x / j;
        sum += term;
        assertEquals(0.019943317853721033, term, 2E-16);
        j = 4;
        term *= x * x / j;
        sum += term;
        assertEquals(0.044872465170872324, term, 2E-16);
        j = 6;
        term *= x * x / j;
        sum += term;
        assertEquals(0.06730869775630849, term, 2E-16);
        j = 8;
        term *= x * x / j;
        sum += term;
        assertEquals(0.07572228497584704, term, 2E-16);
        j = 10;
        term *= x * x / j;
        sum += term;
        assertEquals(0.06815005647826235, term, 2E-16);
        assertEquals(0.28042867064694926, sum, 2E-16);
        double sum2 = 0;
        for (j = 12; j < n; j += 2) {
            term *= x * x / j;
            sum2 += term;
        }
        assertEquals(0.1181924062124705, sum2, 2E-16);
        assertEquals(0.3986210768594198, sum + sum2, 2E-16);

        assertEquals(0.3986210768594198, Momentum.residual(n, x), 2E-16);
    }

    @Test
    public void TestResidualDump() {
        assertTrue(Momentum.residual(Momentum.maxN, (double) Momentum.maxX / Momentum.dividX) < 1E-3);

        System.out.println(System.getProperty("user.dir"));
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/MomentumResidual.txt")) {
            fw.write("x\t");
            for (double x = 0; x <= 16; x += 0.1) {
                fw.write(String.format("%f\t", x));
            }
            fw.write("\n");
            for (int i = 1; i <= 1000; ++i) {
                fw.write(String.format("%d\t", i));
                for (double x = 0; x < 16; x += 0.1) {
                    fw.write(String.format("%e\t", Momentum.residual(i, x)));
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }

    @Test
    public void TestFactorDump() {
        assertTrue(Momentum.residual(Momentum.maxN, (double) Momentum.maxX / Momentum.dividX) < 1E-3);

        System.out.println(System.getProperty("user.dir"));
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/MomentumFactor.txt")) {
            fw.write("x\t");
            for (double x = 0; x <= 16; x += 0.1) {
                fw.write(String.format("%f\t", x));
            }
            fw.write("\n");
            for (int i = 1; i <= 100; ++i) {
                fw.write(String.format("%d\t", i));
                for (double x = 0; x < 16; x += 0.1) {
                    fw.write(String.format("%e\t", Momentum.factor(i, x)));
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }
}
