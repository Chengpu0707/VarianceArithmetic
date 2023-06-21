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
        assertTrue(dbl.toExp(-10));
        assertEquals(-10, dbl.exp);
    }

    @Test
    public void testMaxLimit() {
        alloc(Long.MAX_VALUE / 2);
        assertFalse(dbl.toExp(-1));
        assertEquals(0, dbl.exp);
        assertEquals(Long.MAX_VALUE / 2, dbl.val);

        alloc(Dbl.VAL_EXTRA);
        assertFalse(dbl.toExp(-1));
        assertEquals(0, dbl.exp);
        assertEquals(Dbl.VAL_EXTRA, dbl.val);
    }
    
    @Test
    public void testDirect() {
        alloc(1);
        assertTrue(dbl.toExp(-10));
        assertEquals(-10, dbl.exp);
        assertEquals(1L << 10, dbl.val);

        alloc(8);
        assertTrue(dbl.toExp(-10));
        assertEquals(-10, dbl.exp);
        assertEquals(8L << 10, dbl.val);
    }

    @Test
    public void testIndirect() {
        alloc(Long.MAX_VALUE / 4);
        assertTrue(dbl.toExp(-1));
        assertEquals(-1, dbl.exp);
        assertEquals((Long.MAX_VALUE / 4) * 2, dbl.val);

        alloc(1L << 31);
        assertTrue(dbl.toExp(-30));
        assertEquals(-30, dbl.exp);
        assertEquals(1L << 61, dbl.val);
    }

    @Test
    public void testLimited() {
        alloc(Long.MAX_VALUE / 4);
        assertFalse(dbl.toExp(-2));
        assertEquals(-1, dbl.exp);
        assertEquals((Long.MAX_VALUE / 4) * 2, dbl.val);

        alloc(1L << 31);
        assertFalse(dbl.toExp(-31));
        assertEquals(-30, dbl.exp);
        assertEquals(1L << 61, dbl.val);
    }
}
