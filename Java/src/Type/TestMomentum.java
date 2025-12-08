
package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;



public class TestMomentum {
 
    @Test
    public void OutputNormal() {
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

    public void NormalValues() {
        NormalMomentum norm = new NormalMomentum(5.0);
        assertEquals(norm.bounding, 5.0, 1e-10);
        assertEquals(norm.maxOrder, 448);
        assertEquals(norm.leakage, 2.8665E-7, 1e-10);
        assertEquals(norm.get(0),    1*(1 - 5.733031e-07), 5e-14);
        assertEquals(norm.get(2),    1*(1 - 1.544050e-05), 2e-12);
        assertEquals(norm.get(4),    3*(1 - 1.393338e-04), 3e-11);
        assertEquals(norm.get(6),   15*(1 - 7.588003e-04), 7e-10);
        assertEquals(norm.get(8),  105*(1 - 2.9711805e-03), 2e-9);
        assertEquals(norm.get(10), 945*(1 - 9.1166811e-03), 3e-1);
    }
}
