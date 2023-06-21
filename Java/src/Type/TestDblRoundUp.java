package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.ValueException;

/*
 * Verify Dbl.toExp() is equivalent to mutiple Dbl.onceUp().
 * 
 * Show val+ is equvilent to (val+1)- in rounding.
 * 
 * Some other results on the number of rounding until val == 0:
 *      0- = 0+
 *      1- < 1+
 *      2- < 2+
 *      3- = 3+
 *      4- < 4+
 *      5+ = 5-
 */
public class TestDblRoundUp {
    Dbl dbl;

    private void alloc(long val, boolean rndErr) {
        try {
            dbl = new Dbl(0, false, val, rndErr);
        } catch (ValueException e) {
            fail(e.getMessage());
        }
    }

    private void verify( long val, boolean rndErr, int exp) {
        assertEquals(val, dbl.val);
        assertEquals(rndErr, dbl.rndErr);
        assertEquals(exp, dbl.exp);
    }

    @Test
    public void test0False() {
        alloc(0, false);
        dbl.upOnce();
        verify(0, false, 1);

        alloc(0, false);
        dbl.upBy(1);
        verify(0, false, 1);
    }

    @Test
    public void test0True() {
        alloc(0, true);
        dbl.upOnce();
        verify(0, true, 1);

        alloc(0, true);
        dbl.upBy(1);
        verify(0, true, 1);
    }

    @Test
    public void test1False() {
        alloc(1, false);
        dbl.upOnce();
        verify(0, true, 1);

        alloc(1, false);
        dbl.upBy(1);
        verify(0, true, 1);
    }

    @Test
    public void test1True() {
        alloc(1, true);
        dbl.upOnce();
        verify(1, false, 1);
        dbl.upOnce();
        verify(0, true, 2);

        alloc(1, true);
        dbl.upBy(1);
        verify(1, false, 1);

        alloc(1, true);
        dbl.upBy(2);
        verify(0, true, 2);
    }

    @Test
    public void test2False() {
        alloc(2, false);
        dbl.upOnce();
        verify(1, false, 1);
        dbl.upOnce();
        verify(0, true, 2);

        alloc(2, false);
        dbl.upBy(1);
        verify(1, false, 1);

        alloc(2, false);
        dbl.upBy(2);
        verify(0, true, 2);
    }

    @Test
    public void test2True() {
        alloc(2, true);
        dbl.upOnce();
        verify(1, true, 1);
        dbl.upOnce();
        verify(1, false, 2);
        dbl.upOnce();
        verify(0, true, 3);

        alloc(2, true);
        dbl.upBy(1);
        verify(1, true, 1);

        alloc(2, true);
        dbl.upBy(2);
        verify(1, false, 2);

        alloc(2, true);
        dbl.upBy(3);
        verify(0, true, 3);
    }

    @Test
    public void test3False() {
        alloc(3, false);
        dbl.upOnce();
        verify(1, true, 1);
        dbl.upOnce();
        verify(1, false, 2);
        dbl.upOnce();
        verify(0, true, 3);

        alloc(3, false);
        dbl.upBy(1);
        verify(1, true, 1);

        alloc(3, false);
        dbl.upBy(2);
        verify(1, false, 2);

        alloc(3, false);
        dbl.upBy(3);
        verify(0, true, 3);
    }

    @Test
    public void test3Ture() {
        alloc(3, true);
        dbl.upOnce();
        verify(2, false, 1);
        dbl.upOnce();
        verify(1, false, 2);
        dbl.upOnce();
        verify(0, true, 3);

        alloc(3, true);
        dbl.upBy(1);
        verify(2, false, 1);

        alloc(3, true);
        dbl.upBy(2);
        verify(1, false, 2);

        alloc(3, true);
        dbl.upBy(3);
        verify(0, true, 3);
    }

    @Test
    public void test4False() {
        alloc(4, false);
        dbl.upOnce();
        verify(2, false, 1);
        dbl.upOnce();
        verify(1, false, 2);
        dbl.upOnce();
        verify(0, true, 3);

        alloc(4, false);
        dbl.upBy(1);
        verify(2, false, 1);

        alloc(4, false);
        dbl.upBy(2);
        verify(1, false, 2);

        alloc(4, false);
        dbl.upBy(3);
        verify(0, true, 3);
    }

    @Test
    public void test4True() {
        alloc(4, true);
        dbl.upOnce();
        verify(2, true, 1);
        dbl.upOnce();
        verify(1, true, 2);
        dbl.upOnce();
        verify(1, false, 3);
        dbl.upOnce();
        verify(0, true, 4);

        alloc(4, true);
        dbl.upBy(1);
        verify(2, true, 1);

        alloc(4, true);
        dbl.upBy(2);
        verify(1, true, 2);

        alloc(4, true);
        dbl.upBy(3);
        verify(1, false, 3);

        alloc(4, true);
        dbl.upBy(4);
        verify(0, true, 4);
    }

    @Test
    public void test5False() {
        alloc(5, false);
        dbl.upOnce();
        verify(2, true, 1);
        dbl.upOnce();
        verify(1, true, 2);
        dbl.upOnce();
        verify(1, false, 3);
        dbl.upOnce();
        verify(0, true, 4);

        alloc(5, false);
        dbl.upBy(1);
        verify(2, true, 1);

        alloc(5, false);
        dbl.upBy(2);
        verify(1, true, 2);

        alloc(5, false);
        dbl.upBy(3);
        verify(1, false, 3);

        alloc(5, false);
        dbl.upBy(4);
        verify(0, true, 4);
    }

    @Test
    public void test5True() {
        alloc(5, true);
        dbl.upOnce();
        verify(3, false, 1);
        dbl.upOnce();
        verify(1, true, 2);
        dbl.upOnce();
        verify(1, false, 3);
        dbl.upOnce();
        verify(0, true, 4);

        alloc(5, true);
        dbl.upBy(1);
        verify(3, false, 1);

        alloc(5, true);
        dbl.upBy(2);
        verify(1, true, 2);

        alloc(5, true);
        dbl.upBy(3);
        verify(1, false, 3);

        alloc(5, true);
        dbl.upBy(4);
        verify(0, true, 4);
    }


    @Test
    public void test6False() {
        alloc(6, false);
        dbl.upOnce();
        verify(3, false, 1);
        dbl.upOnce();
        verify(1, true, 2);
        dbl.upOnce();
        verify(1, false, 3);
        dbl.upOnce();
        verify(0, true, 4);

        alloc(6, false);
        dbl.upBy(1);
        verify(3, false, 1);

        alloc(6, false);
        dbl.upBy(2);
        verify(1, true, 2);

        alloc(6, false);
        dbl.upBy(3);
        verify(1, false, 3);

        alloc(6, false);
        dbl.upBy(4);
        verify(0, true, 4);
    }
}
