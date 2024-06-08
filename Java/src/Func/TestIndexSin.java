package Func;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class TestIndexSin {

    @Test
    public void test_neg_rem() {
        // when index is near 0
        assertEquals(0, 1 / 4);     //0 % 1
        assertEquals(1, 1 % 4);
        assertEquals(0, -1 / 4);    //0 % -1
        assertEquals(-1, -1 % 4);

        // when index is near +pi/2
        assertEquals(0, 3 / 4);     //0 % 3
        assertEquals(3, 3 % 4);
        assertEquals(1, 5 / 4);     //0 % 3
        assertEquals(1, 5 % 4);

        // when index is near +pi
        assertEquals(1, 7 / 4);     //0 % 1
        assertEquals(3, 7 % 4);
        assertEquals(2, 9 / 4);     //0 % -1
        assertEquals(1, 9 % 4);

        // when index is near -pi/2
        assertEquals(0, -3 / 4);   //0 % -3
        assertEquals(-3, -3 % 4);
        assertEquals(-1, -5 / 4);   //0 % -3
        assertEquals(-1, -5 % 4);

        // when index is near -pi
        assertEquals(-1, -7 / 4);   //0 % -1
        assertEquals(-3, -7 % 4);
        assertEquals(-2, -9 / 4);   //0 % 1
        assertEquals(-1, -9 % 4);
    }


    @Test
    public void test_get_index() {
        IndexSin indexSin = new IndexSin(3);
        assertEquals(3, indexSin.order());

        assertEquals( 0, indexSin.get_index(0));
        assertEquals( 1, indexSin.get_index(1));
        assertEquals( 2, indexSin.get_index(2));
        assertEquals( 3, indexSin.get_index(3));
        assertEquals( 4, indexSin.get_index(4));
        assertEquals( 3, indexSin.get_index(5));
        assertEquals( 2, indexSin.get_index(6));
        assertEquals( 1, indexSin.get_index(7));
        assertEquals( 0, indexSin.get_index(8));
        assertEquals(-1, indexSin.get_index(9));
        assertEquals(-2, indexSin.get_index(10));
        assertEquals(-3, indexSin.get_index(11));
        assertEquals(-4, indexSin.get_index(12));
        assertEquals(-3, indexSin.get_index(13));
        assertEquals(-2, indexSin.get_index(14));
        assertEquals(-1, indexSin.get_index(15));
        assertEquals( 0, indexSin.get_index(16));
        assertEquals( 1, indexSin.get_index(17));
        assertEquals( 2, indexSin.get_index(18));

        assertEquals(-1, indexSin.get_index(-1));
        assertEquals(-2, indexSin.get_index(-2));
        assertEquals(-3, indexSin.get_index(-3));
        assertEquals(-4, indexSin.get_index(-4));
        assertEquals(-3, indexSin.get_index(-5));
        assertEquals(-2, indexSin.get_index(-6));
        assertEquals(-1, indexSin.get_index(-7));
        assertEquals( 0, indexSin.get_index(-8));
        assertEquals( 1, indexSin.get_index(-9));
        assertEquals( 2, indexSin.get_index(-10));
        assertEquals( 3, indexSin.get_index(-11));
        assertEquals( 4, indexSin.get_index(-12));
        assertEquals( 3, indexSin.get_index(-13));
        assertEquals( 2, indexSin.get_index(-14));
        assertEquals( 1, indexSin.get_index(-15));
        assertEquals( 0, indexSin.get_index(-16));
        assertEquals(-1, indexSin.get_index(-17));
        assertEquals(-2, indexSin.get_index(-18));
    }
 

    static final double q1 = Math.sin(Math.PI/8);
    static final double q2 = Math.sin(Math.PI/4);
    static final double q3 = Math.sin(Math.PI/8*3);

    static final FFT fft = new FFT(FFT.SinSource.IndexSin);

    @Test
    public void testSine() {
        assertEquals(-q3, fft.sin(-5, 4), Math.ulp(q3));

        assertEquals(-1,  fft.sin(-4, 4), Math.ulp(1.0));
        assertEquals(-q3, fft.sin(-3, 4), Math.ulp(q3));
        assertEquals(-q2, fft.sin(-2, 4), Math.ulp(q2));
        assertEquals(-q1, fft.sin(-1, 4), Math.ulp(q1));

        assertEquals(0, fft.sin(0, 4), Double.MIN_VALUE);
        assertEquals(q1, fft.sin(1, 4), Math.ulp(q1));
        assertEquals(q2, fft.sin(2, 4), Math.ulp(q2));
        assertEquals(q3, fft.sin(3, 4), Math.ulp(q3));

        assertEquals(1, fft.sin(4, 4), Math.ulp(1.0));
        assertEquals(q3, fft.sin(5, 4), Math.ulp(q3));
        assertEquals(q2, fft.sin(6, 4), Math.ulp(q2));
        assertEquals(q1, fft.sin(7, 4), Math.ulp(q1));

        assertEquals(0, fft.sin(8, 4), Double.MIN_VALUE);
        assertEquals(-q1, fft.sin(9, 4), Math.ulp(q1));
        assertEquals(-q2, fft.sin(10, 4), Math.ulp(q2));
        assertEquals(-q3, fft.sin(11, 4), Math.ulp(q3));
        assertEquals(-1, fft.sin(12, 4), Math.ulp(1.0));

        assertEquals(-q3, fft.sin(13, 4), Math.ulp(q3));
        assertEquals(-q2, fft.sin(14, 4), Math.ulp(q2));
        assertEquals(-q1, fft.sin(15, 4), Math.ulp(q1));
        assertEquals(0, fft.sin(16, 4), Double.MIN_VALUE);

        assertEquals(q1, fft.sin(17, 4), Math.ulp(q1));
        assertEquals(q2, fft.sin(18, 4), Math.ulp(q2));
        assertEquals(q3, fft.sin(19, 4), Math.ulp(q3));
        
        assertEquals(1, fft.sin(20, 4), Math.ulp(1.0));
    }

    @Test
    public void testLargeIndexSin() {
        assertEquals(fft.sin(16, IndexSin.MAX_ORDER), fft.sin(268451838, 15), Math.ulp(1.0));
        assertEquals(0, fft.sin(-2147483648L, 17), Math.ulp(1.0));
        assertEquals(fft.sin(49152L, 18), fft.sin(2147532800L, 18), Math.ulp(1.0));
    }

    @Test
    public void testCosine() {
        assertEquals(-q1, fft.cos(-5, 4), Math.ulp(q1));

        assertEquals(0, fft.cos(-4, 4), Double.MIN_VALUE);
        assertEquals(q1, fft.cos(-3, 4), Math.ulp(q1));
        assertEquals(q2, fft.cos(-2, 4), Math.ulp(q2));
        assertEquals(q3, fft.cos(-1, 4), Math.ulp(q3));

        assertEquals(1, fft.cos(0, 4), Math.ulp(1.0));
        assertEquals(q3, fft.cos(1, 4), Math.ulp(q3));
        assertEquals(q2, fft.cos(2, 4), Math.ulp(q2));
        assertEquals(q1, fft.cos(3, 4), Math.ulp(q1));

        assertEquals(0, fft.cos(4, 4), Double.MIN_VALUE);
        assertEquals(-q1, fft.cos(5, 4), Math.ulp(q1));
        assertEquals(-q2, fft.cos(6, 4), Math.ulp(q2));
        assertEquals(-q3, fft.cos(7, 4), Math.ulp(q3));

        assertEquals(-1, fft.cos(8, 4), Math.ulp(1.0));
        assertEquals(-q3, fft.cos(9, 4), Math.ulp(q3));
        assertEquals(-q2, fft.cos(10, 4), Math.ulp(q2));
        assertEquals(-q1, fft.cos(11, 4), Math.ulp(q1));

        assertEquals(0, fft.cos(12, 4), Double.MIN_VALUE);
        assertEquals(q1, fft.cos(13, 4), Math.ulp(q1));
        assertEquals(q2, fft.cos(14, 4), Math.ulp(q2));
        assertEquals(q3, fft.cos(15, 4), Math.ulp(q3));

        assertEquals(1, fft.cos(16, 4), Math.ulp(1.0));
        assertEquals(q3, fft.cos(17, 4), Math.ulp(q3));
        assertEquals(q2, fft.cos(18, 4), Math.ulp(q2));
        assertEquals(q1, fft.cos(19, 4), Math.ulp(q1));

        assertEquals(0, fft.cos(20, 4), Double.MIN_VALUE);
    }

    
}
