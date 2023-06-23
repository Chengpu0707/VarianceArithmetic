package Type;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class TestDblMsb {
    @Test
    public void test() {
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
        assertEquals(62, Dbl.msb(1L << 62));
        assertEquals(62, Dbl.msb(Long.MAX_VALUE));
    }
}