
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.FileWriter;
import java.io.IOException;


public class TestMomentum {

    void validate(VarDbl var, double denom, double valueDelta, double uncertainty) {
        final double val = var.value() / denom;
        assertEquals(1, val, valueDelta);
        assertEquals(uncertainty, var.uncertainty(), Math.ulp(uncertainty));
    }

     @Test
    public void TestFactor() {
        validate(Momentum.factor(0), 1, 0.99999947194737793, 1.0510612137199506E-17);
        validate(Momentum.factor(2), 1, 0.99998569318184616, 1.0096418147444105E-17);
        validate(Momentum.factor(4), 3, 0.99987013368935596, 2.998915588392244E-17);
        validate(Momentum.factor(6), 15, 0.99928863789470035, 1.4926804198403977E-16);
        validate(Momentum.factor(8), 105, 0.99719860134891214, 1.0409495599106208E-15);
        validate(Momentum.factor(10), 945, 0.99135593485973217, 9.33728878727979E-15);
      
        assertEquals(0, Momentum.factor(1).value(), 0);
        assertEquals(0, Momentum.factor(3).value(), 0);
        assertEquals(0, Momentum.factor(5).value(), 0);
        assertEquals(0, Momentum.factor(7).value(), 0);
        assertEquals(0, Momentum.factor(9).value(), 0);
    }

    /*
     * Generate data for IPyNb/Momentum.ipynb.
     */
    @Test
    public void TestFactorDump() {
        System.out.println(System.getProperty("user.dir"));
        try (
            final FileWriter fw = new FileWriter("./Java/Output/Momentum.txt")) {
            fw.write("2n\t(2n-1)!!\tValue\tUncertainty\n");
            double doubleFactor = 1;
            for (int i = 2; i <= Momentum.MAX_FACTOR; i += 2, doubleFactor *= (i - 1)) {
                final VarDbl fac = Momentum.factor(i);
                fw.write(String.format("%d\t%e\t%e\t%e\n", 
                         i, doubleFactor, fac.value(), fac.uncertainty()));
            }
            fw.close();
        } catch (IOException e) {
            fail(e.getMessage());
        } 
    }
}
