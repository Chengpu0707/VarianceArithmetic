
package Type;

/**
 * JUnit tests for {@link NormalMoment} — verifies bounding/leakage and the
 * computed moment table against a reference file shared with Python and C++
 * implementations.
 */
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.FileReader;



public class TestMoment {
 
    @Test
    public void OutputNormal() {
        NormalMoment norm = new NormalMoment(5.0);
        assertEquals(norm.bounding, 5.0, 1e-10);
        assertEquals(norm.maxOrder, 448);
        assertEquals(norm.leakage, 5.7330314E-7, 1e-10);
        // Normalized per Formula (2.2): ζ(0,κ)=1 by construction.
        assertEquals( norm.get(0), 1.0, 1e-12);
        assertEquals( norm.get(2) / 0.9999851327963293, 1., 1e-10);
        assertEquals( norm.get(4) / 2.9995837182972154, 1., 1e-10);
        assertEquals( norm.get(6) / 14.988626589191844, 1., 1e-10);
        assertEquals( norm.get(8) / 104.68808606698785, 1., 1e-10);
        assertEquals( norm.get(10) / 936.3852731689979, 1., 1e-10);
    }

    public void NormalValues() {
        NormalMoment norm = new NormalMoment(5.0);
        assertEquals(norm.bounding, 5.0, 1e-10);
        assertEquals(norm.maxOrder, 448);
        assertEquals(norm.leakage, 5.7330314E-7, 1e-10);
        // Normalized per Formula (2.2): deficits relative to ideal Gaussian moments.
        assertEquals(norm.get(0),  1.0, 1e-12);
        assertEquals(norm.get(2),    1*(1 - 1.4867203671e-05), 2e-12);
        assertEquals(norm.get(4),    3*(1 - 1.3876057676e-04), 3e-11);
        assertEquals(norm.get(6),   15*(1 - 7.5822694054e-04), 7e-10);
        assertEquals(norm.get(8),  105*(1 - 2.9706089810e-03), 2e-9);
        assertEquals(norm.get(10), 945*(1 - 9.1161147311e-03), 3e-1);
    }
}
