
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;


public class TestMomentum {

    @Test
    public void TestDoubleFactorial() {
        assertEquals(1, Momentum.doubleFactorial(0), 0);
        assertEquals(1, Momentum.doubleFactorial(1), 0);
        assertEquals(2, Momentum.doubleFactorial(2), 0);
        assertEquals(3, Momentum.doubleFactorial(3), 0);
        assertEquals(8, Momentum.doubleFactorial(4), 0);
        assertEquals(15, Momentum.doubleFactorial(5), 0);
        assertEquals(2.7253921397507295E78, Momentum.doubleFactorial(99), 1E69);
        assertEquals(1.0898143681335199E154, Momentum.doubleFactorial(170), 1E145);
        assertEquals(1.863582569508319E156, Momentum.doubleFactorial(170)*171, 1E147);
        assertEquals(1.183050330245448E188, Momentum.doubleFactorial(200), 1E179);
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
        assertEquals(  1 * 0.99999947194737793, Momentum.factor(0), Math.ulp(1));
        assertEquals(  1 * 0.99998569318184616, Momentum.factor(2), 2e-6);
        assertEquals(  3 * 0.99987013368935596, Momentum.factor(4), 3e-5);
        assertEquals( 15 * 0.99928863789470035, Momentum.factor(6), 1e-3);
        assertEquals(105 * 0.99719860134891214, Momentum.factor(8), 0.1);
        assertEquals(945 * 0.99135593485973217, Momentum.factor(10), 0.5); 
    
        assertEquals(0, Momentum.factor(1), 0);
        assertEquals(0, Momentum.factor(3), 0);
        assertEquals(0, Momentum.factor(5), 0);
        assertEquals(0, Momentum.factor(7), 0);
        assertEquals(0, Momentum.factor(9), 0);
    }

    @Test
    public void TestFactorDump() {
        System.out.println(System.getProperty("user.dir"));
        try (
            final FileWriter fw = new FileWriter("./Java/Output/MomentumFactor.txt")) {
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
