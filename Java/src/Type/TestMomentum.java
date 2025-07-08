
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;



public class TestMomentum {
 
    @Test
    public void TestNormal() {
        NormalMomentum norm = new NormalMomentum(5.0);
        assertEquals(norm.bounding, 5.0, 1e-10);
        assertEquals(norm.maxOrder, 448);
        assertEquals(norm.leakage, 2.8665E-7, 1e-10);

        System.out.println(System.getProperty("user.dir"));
        try (final BufferedReader fr = new BufferedReader(new FileReader("./Python/NormalMomentum_5.0.txt"))) {
            String line = fr.readLine().strip();
            assertEquals(line, "n\tMomentum\tBounding:\t5.0");
            int i = 0;
            while ((line = fr.readLine()) != null) {
                line = line.strip();
                final String[] sWord = line.split("\t");
                assertEquals(i*2, Integer.parseInt(sWord[0]));
                final double value = Double.parseDouble(sWord[1]);
                assertEquals(norm.get(i*2) /value, 1.0, 1e-10);
                assertEquals(norm.get(i*2 + 1), 0, 1e-17);
                ++i;
            }
        } catch (Throwable e) {
            fail(e.getMessage());
        } 
    }
}
