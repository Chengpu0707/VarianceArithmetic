
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
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
    public void TestFactorNorm() {
        assertEquals(1.0, Momentum.factor(2, 6, false), 1E-6);
        assertEquals(1.0, Momentum.factor(2, 5, false), 2E-4);
        assertEquals(1.0, Momentum.factor(2, 4, false), 2E-3);

        assertEquals(1, Momentum.factor(2, 5, true), 0);
        assertEquals(1, Momentum.factor(4, 5, true)/3, 2E-4);
        assertEquals(1, Momentum.factor(6, 5, true)/15, 3E-3);
        assertEquals(1, Momentum.factor(8, 5, true)/105, 3E-3);
        assertEquals(1, Momentum.factor(10, 5, true)/945, 1E-2);
    }


    @Test
    public void TestFactor() {
        final double s = 1.75;
        final int n = 2;
        assertEquals(0.617911, Momentum.factor(n, s, false), 1E-6);
        final double pdf = 2 * s * s * Momentum.pdf(s);
        assertEquals(0.528449, pdf, 1E-6);

        double term = pdf * s / (n + 1);
        double sum = 0;
        double next = term;
        for (int j = n + 3; sum < next; j += 2) {
            sum = next;
            term *= s * s / j;
            next = sum + term;
        }
        assertEquals(sum, Momentum.factor(n, s, false), 1E-6);
        assertEquals(Momentum.factor(n+2, s, false), (n+1) * sum - pdf*s, 1E-6);

        sum = 0;
        term = pdf * s /(n + 1);
        sum += term;
        assertEquals(0.308262, term, 1E-6);
        term *= s * s /(n + 3);
        assertEquals(0.188810, term, 1E-6);
        sum += term;
        term *= s * s /(n + 5);
        assertEquals(0.082604, term, 1E-6);
        term *= s * s /(n + 7);
        assertEquals(0.028108, term, 1E-6);
        sum += term;
        assertEquals(sum + 0.092731, Momentum.factor(n, s, false), 1E-6);
    }

    @Test
    public void TestFactorDump() {
        System.out.println(System.getProperty("user.dir"));
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/MomentumFactor.txt")) {
            final double maxS = (double) Momentum.maxS / Momentum.dividS;
            fw.write("s\t(2n-1)!!\t");
            for (double s = 0.1; s < maxS; s += 0.1) {
                fw.write(String.format("%f\t", s));
            }
            fw.write("\n");
            for (int i = 2; i <= Momentum.maxN; i += 2) {
                fw.write(String.format("%d\t%e\t", i, Momentum.doubleFactorial(i - 1)));
                for (double s = 0.1; s < maxS; s += 0.1) {
                    fw.write(String.format("%e\t", Momentum.factor(i, s, false)));
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }
}
