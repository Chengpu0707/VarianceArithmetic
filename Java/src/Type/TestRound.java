/*
 * Test Round on:
 *  *) The agreement of upOnce() and upBy()
 */
package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class TestRound {
    Round rnd;
    
    @Test
    public void test0() {
        rnd = new Round(0);
        assertFalse(rnd.upOnce());
        assertEquals(0, rnd.val);
        assertEquals(false, rnd.rndErr);
    }

    @Test
    public void test1() {
        rnd = new Round(1);
        assertTrue(rnd.upOnce());
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(1);
        rnd.upBy(1);
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);
    }

    @Test
    public void test2() {
        rnd = new Round(2);
        assertFalse(rnd.upOnce());
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);
        assertTrue(rnd.upOnce());
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(2);
        assertFalse(rnd.upBy(1));
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);

        rnd = new Round(2);
        assertTrue(rnd.upBy(2));
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);
    }

    @Test
    public void test3() {
        rnd = new Round(3);
        rnd.upOnce();
        assertEquals(1, rnd.val);
        assertEquals(true, rnd.rndErr);
        rnd.upOnce();
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);
        rnd.upOnce();
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(3);
        rnd.upBy(1);
        assertEquals(1, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(3);
        rnd.upBy(2);
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);

        rnd = new Round(3);
        rnd.upBy(3);
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);
    }

    @Test
    public void test4() {
        rnd = new Round(4);
        rnd.upOnce();
        assertEquals(2, rnd.val);
        assertEquals(false, rnd.rndErr);
        rnd.upOnce();
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);
        rnd.upOnce();
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(4);
        rnd.upBy(1);
        assertEquals(2, rnd.val);
        assertEquals(false, rnd.rndErr);

        rnd = new Round(4);
        rnd.upBy(2);
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);

        rnd = new Round(4);
        rnd.upBy(3);
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);
    }

    @Test
    public void test5() {
        rnd = new Round(5);
        rnd.upOnce();
        assertEquals(2, rnd.val);
        assertEquals(true, rnd.rndErr);
        rnd.upOnce();
        assertEquals(1, rnd.val);
        assertEquals(true, rnd.rndErr);
        rnd.upOnce();
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);
        rnd.upOnce();
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(5);
        rnd.upBy(1);
        assertEquals(2, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(5);
        rnd.upBy(2);
        assertEquals(1, rnd.val);
        assertEquals(true, rnd.rndErr);

        rnd = new Round(5);
        rnd.upBy(3);
        assertEquals(1, rnd.val);
        assertEquals(false, rnd.rndErr);

        rnd = new Round(5);
        rnd.upBy(4);
        assertEquals(0, rnd.val);
        assertEquals(true, rnd.rndErr);
    }

    @Test
    public void testDownTo() {
        rnd = new Round(6);
        assertFalse(rnd.upBy(-2));
        assertEquals(6 << 2, rnd.val);
        assertEquals(false, rnd.rndErr);
    }

    @Test
    public void testLargeUp() {
        assertEquals( 4611686018427387904L, 1L << 62);
        assertEquals(-9223372036854775808L, 1L << 63);
        // the mask function no longer works
        assertEquals(1, 1L << 64);
        assertEquals(2, 1L << 65);

        rnd = new Round(Dbl.DOUBLE_VAL_EXTRA);
        assertTrue(rnd.upBy(Long.SIZE));
        assertEquals(0, rnd.val);
        rnd = new Round(Dbl.DOUBLE_VAL_EXTRA);
        assertTrue(rnd.upBy(Long.SIZE + 1));
        assertEquals(0, rnd.val);
    }

    @Test
    public void testLargeDown() {
        assertEquals(1, 1L << Long.SIZE);
        assertEquals(2, 1L << (Long.SIZE + 1));

        rnd = new Round(1);
        assertFalse(rnd.upBy(-Long.SIZE));
        assertEquals(1, rnd.val);
        rnd = new Round(1);
        assertFalse(rnd.upBy(-Long.SIZE-1));
        assertEquals(2, rnd.val);
    }
}
