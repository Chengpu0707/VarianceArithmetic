
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
    public void TestResidual() {
        System.out.println(System.getProperty("user.dir"));
        try (
            final FileWriter fw = new FileWriter("C:/Users/Cheng/Documents/Proj/VarianceArithemtic/Java/Output/MomentumResidual.txt")) {
            fw.write("x\t");
            for (double x = 0; x <= 12; x += 0.1) {
                fw.write(String.format("%f\t", x));
            }
            fw.write("\n");
            for (int i = 1; i <= 50; ++i) {
                fw.write(String.format("%d\t", i));
                for (double x = 0; x < 12; x += 0.1) {
                    fw.write(String.format("%e\t", Momentum.residual(i, x)));
                }
                fw.write("\n");
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }
}
