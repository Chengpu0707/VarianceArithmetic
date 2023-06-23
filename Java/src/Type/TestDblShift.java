package Type;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import Type.IReal.ValueException;

public class TestDblShift {
    Dbl dbl;
    Dbl res;

    void test(int exp, long val, int shift) {
        try {
            dbl = new Dbl(exp, false, val, false);
            res = new Dbl(dbl);
            assertTrue(res.shift(shift));
            if (shift >= 0) {
                assertEquals(dbl.toDouble() * (1L << shift), res.toDouble(), 0);
            } else {
                assertEquals(dbl.toDouble() / (1L << -shift), res.toDouble(), 0);
            }
        } catch (ValueException e) {
            fail(e.getMessage());
        }
    }

    void testLimit(int exp, long val, int shift) {
        try {
            dbl = new Dbl(exp, false, val, false);
            res = new Dbl(dbl);
            assertFalse(res.shift(shift));
            if (shift >= 0) {
                assertEquals(Integer.MAX_VALUE, res.exp());
            } else {
                assertEquals(Integer.MIN_VALUE, res.exp());
            }
        } catch (ValueException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void test0() {
        test(0, 0, 0);
        test(0, 0, 10);
        test(0, 0, -10);
        test(10, 0, 0);
        test(10, 0, 10);
        test(10, 0, -10);
        test(-10, 0, 0);
        test(-10, 0, 10);
        test(-10, 0, -10);
        testLimit(Integer.MAX_VALUE - 5, 0, 10);
        testLimit(Integer.MIN_VALUE + 5, 0, -10);
    }

    @Test
    public void test1() {
        test(0, 1, 0);
        test(0, 1, 10);
        test(0, 1, -10);
        test(10, 1, 0);
        test(10, 1, 10);
        test(10, 1, -10);
        test(-10, 1, 0);
        test(-10, 1, 10);
        test(-10, 1, -10);
        testLimit(Integer.MAX_VALUE - 5, 0, 10);
        testLimit(Integer.MIN_VALUE + 5, 0, -10);
    }

}
