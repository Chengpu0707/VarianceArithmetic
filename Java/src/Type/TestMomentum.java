
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
        assertEquals( norm.get(0) / 0.999999426696856, 1., 1e-10);
        assertEquals( norm.get(2) / 0.999984559501709, 1., 1e-10);
        assertEquals( norm.get(4) / 2.99958199862644, 1., 1e-10);
        assertEquals( norm.get(6) / 14.9886179961651, 1., 1e-10);
        assertEquals( norm.get(8) / 104.688026048979, 1., 1e-10);
        assertEquals( norm.get(10) / 936.384736336377, 1., 1e-10);
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
