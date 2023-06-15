package Type;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/*
 * Some properties of double type
 */
public class TestDouble {
    @Test
    public void testNaN() {
        assertTrue(Double.isNaN(Double.NaN));
        assertTrue(Double.isNaN(-Double.NaN));
        assertFalse(Double.isInfinite(Double.NaN));
        assertFalse(Double.isInfinite(-Double.NaN));
        assertEquals(Double.NaN, Double.NaN, 0);
    }

    @Test
    public void testInfinity() {
        assertTrue(Double.isInfinite(Double.POSITIVE_INFINITY));
        assertTrue(Double.isInfinite(Double.NEGATIVE_INFINITY));
        assertFalse(Double.isNaN(Double.POSITIVE_INFINITY));
        assertFalse(Double.isNaN(Double.NEGATIVE_INFINITY));

        assertEquals(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 0);
        assertEquals(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 0);
        assertNotEquals(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 0);
        assertNotEquals(Double.POSITIVE_INFINITY, Double.NaN, 0);
        assertNotEquals(Double.NEGATIVE_INFINITY, Double.NaN, 0);

        assertTrue(+Double.MAX_VALUE < Double.POSITIVE_INFINITY);
        assertTrue(Double.NEGATIVE_INFINITY < -Double.MAX_VALUE);
    }

    @Test
    public void testMinMax() {
        assertEquals(Double.MIN_NORMAL, 4 / Double.MAX_VALUE, Double.MIN_NORMAL * 2E-16);
        assertEquals(Double.MAX_VALUE / 4, 1 / Double.MIN_NORMAL, Double.MAX_VALUE / 4 * 2E-16);
    }

    @Test
    public void testOverflow() {
        assertEquals(+Double.MAX_VALUE, +Double.MAX_VALUE + Double.MAX_VALUE * 0.5E-16, 0);
        assertEquals(-Double.MAX_VALUE, -Double.MAX_VALUE - Double.MAX_VALUE * 0.5E-16, 0);
        assertEquals(Double.POSITIVE_INFINITY, +Double.MAX_VALUE + Double.MAX_VALUE * 1E-16, 0);
        assertEquals(Double.NEGATIVE_INFINITY, -Double.MAX_VALUE - Double.MAX_VALUE * 1E-16, 0);

        assertEquals(+Double.MAX_VALUE, +Double.MAX_VALUE * (1 + 1E-16), 0);
        assertEquals(-Double.MAX_VALUE, -Double.MAX_VALUE - (1 + 1E-16), 0);
        assertEquals(Double.POSITIVE_INFINITY, +Double.MAX_VALUE * (1 + 2E-16), 0);
        assertEquals(Double.NEGATIVE_INFINITY, -Double.MAX_VALUE * (1 + 2E-16), 0);

        // the long overflow 
        assertEquals(Long.MIN_VALUE, Long.MAX_VALUE + 1L);
        assertEquals(Long.MAX_VALUE, Long.MIN_VALUE - 1L);
        assertEquals(-2L, Long.MAX_VALUE * 2L);
        assertEquals( 0L, Long.MIN_VALUE * 2L);

        // the int overflow 
        assertEquals(Integer.MIN_VALUE, Integer.MAX_VALUE + 1);
        assertEquals(Integer.MAX_VALUE, Integer.MIN_VALUE - 1);
        assertEquals(-2, Integer.MAX_VALUE * 2);
        assertEquals( 0L, Integer.MIN_VALUE * 2);
    }

    @Test
    public void testDividZero() {
        assertEquals(Double.POSITIVE_INFINITY, +1.0 / 0, 0);
        assertEquals(Double.NEGATIVE_INFINITY, -1.0 / 0, 0);

        // int divid zero will cause an exception
        try {
            int i = 1 / 0;
            fail();
        } catch( Exception e ) {
        }
        try {
            long i = 1L / 0;
            fail();
        } catch( Exception e ) {
        }
    }

    @Test
    public void testPowExp0() {
        assertEquals(1, Math.pow( 2, 0 ), 0);
        assertEquals(1, Math.pow( -2, 0 ), 0);
        assertEquals(1, Math.pow( Double.NaN, 0 ), 0);
        assertEquals(1, Math.pow( Double.POSITIVE_INFINITY, 0 ), 0);
        assertEquals(1, Math.pow( Double.NEGATIVE_INFINITY, 0 ), 0);
    }

    @Test
    public void testPowExp1() {
        assertEquals(2, Math.pow( 2, 1 ), 0);
        assertEquals(-2, Math.pow( -2, 1 ), 1);
        assertEquals(Double.NaN, Math.pow( Double.NaN, 1 ), 0);
        assertEquals(Double.POSITIVE_INFINITY, Math.pow( Double.POSITIVE_INFINITY, 1 ), 0);
        assertEquals(Double.NEGATIVE_INFINITY, Math.pow( Double.NEGATIVE_INFINITY, 1 ), 0);
    }

    @Test
    public void testSqrt() {
        assertEquals(0, Math.sqrt( 0 ), 0);
        assertEquals(2, Math.sqrt( 4 ), 0);
        assertEquals(Double.NaN, Math.sqrt( -4 ), 0);
        assertEquals(Double.NaN, Math.sqrt( Double.NaN ), 0);
        assertEquals(Double.POSITIVE_INFINITY, Math.sqrt( Double.POSITIVE_INFINITY ), 0);
        assertEquals(Double.NaN, Math.sqrt( Double.NEGATIVE_INFINITY ), 0);

        assertEquals(0, Math.pow( 0, 0.5 ), 0);
        assertEquals(2, Math.pow( 4, 0.5 ), 0);
        assertEquals(Double.NaN, Math.pow( -4, 0.5 ), 0);
        assertEquals(Double.NaN, Math.pow( Double.NaN, 0.5 ), 0);
        assertEquals(Double.POSITIVE_INFINITY, Math.pow( Double.POSITIVE_INFINITY, 0.5 ), 0);
        assertEquals(Double.POSITIVE_INFINITY, Math.pow( Double.NEGATIVE_INFINITY, 0.5 ), 0);
            // Math.sqrt( Double.NEGATIVE_INFINITY ) and Math.pow( Double.NEGATIVE_INFINITY, 0.5 ) disagree
    }
    
    @Test
    public void testInvSqrt() {
        assertEquals(Double.POSITIVE_INFINITY, Math.pow( 0, -0.5 ), 0);
        assertEquals(0.5, Math.pow( 4, -0.5 ), 0);
        assertEquals(Double.NaN, Math.pow( -4, -0.5 ), 0);
        assertEquals(Double.NaN, Math.pow( Double.NaN, -0.5 ), 0);
    }

    @Test
    public void testWhenInfinityBecomeZero() {
        assertEquals(0, +1.0 / Double.POSITIVE_INFINITY, 0);
        assertEquals(0, +1.0 / Double.NEGATIVE_INFINITY, 0);
        assertEquals(0, Math.pow( Double.POSITIVE_INFINITY, -0.5 ), 0);
        assertEquals(0, Math.pow( Double.NEGATIVE_INFINITY, -0.5 ), 0);
    }
    
    @Test
    public void testRoundingError() {
        assertEquals(1, Math.sqrt(2) / Math.sqrt(2), 1E-16);
        assertEquals(1, Math.sqrt(2) * Math.sqrt(2) / 2, 3E-16);
        assertEquals(1, Math.sqrt(2) / Math.sqrt(0.5) / 2, 2E-16);
        assertEquals(1, Math.sqrt(0.5) * Math.sqrt(0.5) / 0.5, 3E-16);

        assertEquals(1, Math.sqrt(5) * Math.sqrt(5) / 5, 3E-16);

        assertEquals(1, 2.0/3 + 1.0/3, 1E-16);
        assertEquals(1.0/3, 2.0/3 - 1.0/3, 1E-16);
        assertEquals(1, 1.0/11 + 10.0/11, 1E-16);
    }
};
  
  

