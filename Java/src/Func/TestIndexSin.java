package Func;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class TestIndexSin {
    private IndexSin indexSin = new IndexSin(SinSource.Quart);

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
        assertEquals( 0, indexSin.get_index(0, 3));
        assertEquals( 1, indexSin.get_index(1, 3));
        assertEquals( 2, indexSin.get_index(2, 3));
        assertEquals( 3, indexSin.get_index(3, 3));
        assertEquals( 4, indexSin.get_index(4, 3));
        assertEquals( 3, indexSin.get_index(5, 3));
        assertEquals( 2, indexSin.get_index(6, 3));
        assertEquals( 1, indexSin.get_index(7, 3));
        assertEquals( 0, indexSin.get_index(8, 3));
        assertEquals(-1, indexSin.get_index(9, 3));
        assertEquals(-2, indexSin.get_index(10, 3));
        assertEquals(-3, indexSin.get_index(11, 3));
        assertEquals(-4, indexSin.get_index(12, 3));
        assertEquals(-3, indexSin.get_index(13, 3));
        assertEquals(-2, indexSin.get_index(14, 3));
        assertEquals(-1, indexSin.get_index(15, 3));
        assertEquals( 0, indexSin.get_index(16, 3));
        assertEquals( 1, indexSin.get_index(17, 3));
        assertEquals( 2, indexSin.get_index(18, 3));

        assertEquals(-1, indexSin.get_index(-1, 3));
        assertEquals(-2, indexSin.get_index(-2, 3));
        assertEquals(-3, indexSin.get_index(-3, 3));
        assertEquals(-4, indexSin.get_index(-4, 3));
        assertEquals(-3, indexSin.get_index(-5, 3));
        assertEquals(-2, indexSin.get_index(-6, 3));
        assertEquals(-1, indexSin.get_index(-7, 3));
        assertEquals( 0, indexSin.get_index(-8, 3));
        assertEquals( 1, indexSin.get_index(-9, 3));
        assertEquals( 2, indexSin.get_index(-10, 3));
        assertEquals( 3, indexSin.get_index(-11, 3));
        assertEquals( 4, indexSin.get_index(-12, 3));
        assertEquals( 3, indexSin.get_index(-13, 3));
        assertEquals( 2, indexSin.get_index(-14, 3));
        assertEquals( 1, indexSin.get_index(-15, 3));
        assertEquals( 0, indexSin.get_index(-16, 3));
        assertEquals(-1, indexSin.get_index(-17, 3));
        assertEquals(-2, indexSin.get_index(-18, 3));
    }
 

    static final double q1 = Math.sin(Math.PI/8);
    static final double q2 = Math.sin(Math.PI/4);
    static final double q3 = Math.sin(Math.PI/8*3);

    private void testSine(IndexSin index, double delta) {
        assertEquals(0, index.sin(0, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);
        assertEquals(q1, index.sin(1, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(q2, index.sin(2, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q3, index.sin(3, 3).value(), (delta > 0)? delta : Math.ulp(q3));

        assertEquals(1, index.sin(4, 3).value(), (delta > 0)? delta : Math.ulp(1.0));
        assertEquals(q3, index.sin(5, 3).value(), (delta > 0)? delta : Math.ulp(q3));
        assertEquals(q2, index.sin(6, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q1, index.sin(7, 3).value(), (delta > 0)? delta : Math.ulp(q1));

        assertEquals(0, index.sin(8, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);
        assertEquals(-q1, index.sin(9, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(-q2, index.sin(10, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(-q3, index.sin(11, 3).value(),(delta > 0)? delta :  Math.ulp(q3));
        assertEquals(-1, index.sin(12, 3).value(), (delta > 0)? delta : Math.ulp(1.0));

        assertEquals(-q3, index.sin(13, 3).value(), (delta > 0)? delta : Math.ulp(q3));
        assertEquals(-q2, index.sin(14, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(-q1, index.sin(15, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(0, index.sin(16, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);

        assertEquals(q1, index.sin(17, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(q2, index.sin(18, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q3, index.sin(19, 3).value(), (delta > 0)? delta : Math.ulp(q3));       
        assertEquals(1, index.sin(20, 3).value(), (delta > 0)? delta : Math.ulp(1.0));

        assertEquals(-q3, index.sin(-5, 3).value(), (delta > 0)? delta : Math.ulp(q3));

        assertEquals(-1,  index.sin(-4, 3).value(), (delta > 0)? delta : Math.ulp(1.0));
        assertEquals(-q3, index.sin(-3, 3).value(), (delta > 0)? delta : Math.ulp(q3));
        assertEquals(-q2, index.sin(-2, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(-q1, index.sin(-1, 3).value(), (delta > 0)? delta : Math.ulp(q1));
    }

    private void testSine(IndexSin index) {
        testSine(index, 0);
    }

    private void testCosine(IndexSin index, double delta) {
        assertEquals(-q1, index.cos(-5, 3).value(), (delta > 0)? delta : Math.ulp(q1));

        assertEquals(0, index.cos(-4, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);
        assertEquals(q1, index.cos(-3, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(q2, index.cos(-2, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q3, index.cos(-1, 3).value(), (delta > 0)? delta : Math.ulp(q3));

        assertEquals(1, index.cos(0, 3).value(), (delta > 0)? delta : Math.ulp(1.0));
        assertEquals(q3, index.cos(1, 3).value(), (delta > 0)? delta : Math.ulp(q3));
        assertEquals(q2, index.cos(2, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q1, index.cos(3, 3).value(), (delta > 0)? delta : Math.ulp(q1));

        assertEquals(0, index.cos(4, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);
        assertEquals(-q1, index.cos(5, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(-q2, index.cos(6, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(-q3, index.cos(7, 3).value(), (delta > 0)? delta : Math.ulp(q3));

        assertEquals(-1, index.cos(8, 3).value(), (delta > 0)? delta : Math.ulp(1.0));
        assertEquals(-q3, index.cos(9, 3).value(), (delta > 0)? delta : Math.ulp(q3));
        assertEquals(-q2, index.cos(10, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(-q1, index.cos(11, 3).value(), (delta > 0)? delta : Math.ulp(q1));

        assertEquals(0, index.cos(12, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);
        assertEquals(q1, index.cos(13, 3).value(), (delta > 0)? delta : Math.ulp(q1));
        assertEquals(q2, index.cos(14, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q3, index.cos(15, 3).value(), (delta > 0)? delta : Math.ulp(q3));

        assertEquals(1, index.cos(16, 3).value(), (delta > 0)? delta : Math.ulp(1.0));
        assertEquals(q3, index.cos(17, 3).value(), (delta > 0)? delta : Math.ulp(q3));
        assertEquals(q2, index.cos(18, 3).value(), (delta > 0)? delta : Math.ulp(q2));
        assertEquals(q1, index.cos(19, 3).value(), (delta > 0)? delta : Math.ulp(q1));

        assertEquals(0, index.cos(20, 3).value(), (delta > 0)? delta : Double.MIN_VALUE);
    }

    private void testCosine(IndexSin index) {
        testCosine(index, 0);
    }


    @Test
    public void test_Quart() {
        final IndexSin index = new IndexSin(SinSource.Quart);
        testSine(index);
        testCosine(index);
        indexSin.dump("./Java/Output/IndexSin_Quart_18.txt", 18);
    }

    @Test
    public void test_Prec() {
        final IndexSin index = new IndexSin(SinSource.Prec);
        testSine(index);
        testCosine(index);
    }

    @Test
    public void test_Lib() {
        final IndexSin index = new IndexSin(SinSource.Lib);
        testSine(index, Math.ulp(1));
        testCosine(index, Math.ulp(1));
    }

}
