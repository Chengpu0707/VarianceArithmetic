
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;


public class TestMomentum {

    @Test
    public void TestFactor() {
        assertEquals(Momentum.get(0) / 1, 1, 6e-7);
        assertEquals(Momentum.get(2) / 1, 1, 2e-4);
        assertEquals(Momentum.get(4) / 3, 1, 2e-3);
        assertEquals(Momentum.get(6) / 15, 1, 7e-3);
        assertEquals(Momentum.get(8) / 105, 1, 4e-2);
        assertEquals(Momentum.get(10) / 945, 1, 1e-2);
      
        assertEquals(0, Momentum.get(1), 0);
        assertEquals(0, Momentum.get(3), 0);
        assertEquals(0, Momentum.get(5), 0);
        assertEquals(0, Momentum.get(7), 0);
        assertEquals(0, Momentum.get(9), 0);
    }

    /*
     * compare with Python calculation.
     */
    @Test
    public void TestCompare() {
        System.out.println(System.getProperty("user.dir"));
        try (   final BufferedReader fr = new BufferedReader(new FileReader("./Python/Output/NormalMomentum_5.txt"));
                final FileWriter fw = new FileWriter("./Java/Output/NormalMomentum_5.txt")) {
            String line = fr.readLine().strip();
            assertEquals(line, "n\tMomentum\t!!Diff\tSigma=5.0");
            fw.write("2n\tPython\tJava\tError\n");
            for (int i = 0; i < Momentum.MAX_FACTOR; i += 2) {
                line = fr.readLine().strip();
                final String[] sWord = line.strip().split("\t");
                assertEquals(i, Integer.parseInt(sWord[0]));
                final double value = Double.parseDouble(sWord[1]);
                fw.write(String.format("%d\t%e\t%e\t%e\n", i, value, Momentum.get(i), Momentum.get(i)/value - 1));
            }
        } catch (Throwable e) {
            fail(e.getMessage());
        } 
    }
}
