package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.ValueException;


public class TestDblRoundDown {
    Dbl dbl;

    private void alloc(long val) {
        try {
            dbl = new Dbl(0, false, val, true);
        } catch (ValueException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void test0() {
        alloc(0);
        assertTrue(dbl.upBy(-10));
        assertEquals(-10, dbl.exp);
    }

    @Test
    public void testMaxLimit() {
        alloc(Long.MAX_VALUE / 2);
        assertFalse(dbl.upBy(-1, 62));
        assertEquals(0, dbl.exp);
        assertEquals(Long.MAX_VALUE / 2, dbl.val);

        alloc(1L << 62);
        assertFalse(dbl.upBy(-1, 62));
        assertEquals(0, dbl.exp);
        assertEquals(1L << 62, dbl.val);
    }
    
    @Test
    public void testDirect() {
        alloc(1);
        assertTrue(dbl.upBy(-10));
        assertEquals(-10, dbl.exp);
        assertEquals(1L << 10, dbl.val);

        alloc(8);
        assertTrue(dbl.upBy(-10));
        assertEquals(-10, dbl.exp);
        assertEquals(8L << 10, dbl.val);
    }

    @Test
    public void testIndirect() {
        alloc(Long.MAX_VALUE / 4);
        assertTrue(dbl.upBy(-1));
        assertEquals(-1, dbl.exp);
        assertEquals((Long.MAX_VALUE / 4) * 2, dbl.val);

        alloc(1L << 31);
        assertTrue(dbl.upBy(-30));
        assertEquals(-30, dbl.exp);
        assertEquals(1L << 61, dbl.val);
    }

    @Test
    public void testLimited() {
        alloc(Long.MAX_VALUE / 4);
        assertFalse(dbl.upBy(-2));
        assertEquals(-1, dbl.exp);
        assertEquals((Long.MAX_VALUE / 4) * 2, dbl.val);

        alloc(1L << 31);
        assertFalse(dbl.upBy(-31));
        assertEquals(-30, dbl.exp);
        assertEquals(1L << 61, dbl.val);
    }

    @Test
    public void testMSB() {
        assertEquals(63, Dbl.msb(-1));
        assertEquals(0, Dbl.msb(0));
        assertEquals(0, Dbl.msb(1));
        assertEquals(1, Dbl.msb(2));
        assertEquals(1, Dbl.msb(3));
        assertEquals(2, Dbl.msb(4));
        assertEquals(2, Dbl.msb(5));
        assertEquals(10, Dbl.msb(1L << 10));
        assertEquals(Dbl.DOUBLE_EXP_SHIFT, Dbl.msb(Dbl.DOUBLE_VAL_EXTRA));
        assertEquals(Dbl.DOUBLE_EXP_SHIFT, Dbl.msb(Dbl.DOUBLE_VAL_MAX));
        assertEquals(57, Dbl.msb(1L << 57));
        assertEquals(57, Dbl.msb((1L << 57) + (1L << 10)));
        assertEquals(63, Dbl.msb(Long.MAX_VALUE));
    }
}
