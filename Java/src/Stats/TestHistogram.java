package Stats;

import org.junit.Assert;
import org.junit.Test;

public class TestHistogram {
	private static final double TOLERANCE = 3E-16;

	@Test
	public void testDefault() {
		final Histogram histo = new Histogram();
        Assert.assertEquals(10, histo.divids());
        Assert.assertEquals(3, histo.maxRange(), 1E-16);
        final double[] sExp = new double[61];
		Assert.assertTrue( null == histo.histo() );
		Assert.assertTrue( histo.accum( 0, 1 ) );
        sExp[30] = 1;
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( 0.01, 2 ) );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( -0.01, 3 ) );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		
		Assert.assertFalse( histo.accum( 10, 4 ) );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );

		Assert.assertTrue( histo.accum( -0.5, 5 ) );
		Assert.assertEquals( 0, histo.lower() );
        Assert.assertEquals( 1, histo.upper() );
       sExp[30] = 3.0/4;
        sExp[25] = 1.0/4;
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( 0.295, 6 ) );
        sExp[30] = 3.0/5;
        sExp[25] = 1.0/5;
        sExp[33] = 1.0/5;
		Assert.assertEquals( 0, histo.lower() );
        Assert.assertEquals( 1, histo.upper() );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertArrayEquals( new double[] {0, 0, 0, 3.0/5, 0, 0, 1.0/5}, 
				histo.histo(3), TOLERANCE );
		Assert.assertArrayEquals( new double[] {0, 1.0/5, 0, 0, 0, 0, 3.0/5, 0, 0, 1.0/5, 0, 0, 0}, 
				histo.histo(6), TOLERANCE );
	}
}
