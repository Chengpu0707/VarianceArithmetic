package Stats;

import org.junit.Assert;
import org.junit.Test;

public class TestHistogram {
	private static final double TOLERANCE = 3E-16;

	@Test
	public void testDefault() {
		final Histogram histo = new Histogram();
		Assert.assertTrue( null == histo.histo() );
		Assert.assertTrue( histo.accum( 0 ) );
		Assert.assertArrayEquals( new double[] {0}, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( 0.01 ) );
		Assert.assertArrayEquals( new double[] {1}, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( -0.01 ) );
		Assert.assertArrayEquals( new double[] {1}, histo.histo(), TOLERANCE );
		
		Assert.assertFalse( histo.accum( 10 ) );

		Assert.assertTrue( histo.accum( -0.5 ) );
		Assert.assertEquals( 5, histo.actRange() );
		Assert.assertArrayEquals( new double[] {1.0/4, 0, 0, 0, 0, 3.0/4, 0, 0, 0, 0, 0}, 
				histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( 0.295 ) );
		Assert.assertEquals( 5, histo.actRange() );
		Assert.assertArrayEquals( new double[] {1.0/5, 0, 0, 0, 0, 3.0/5, 0, 0, 1.0/5, 0, 0}, 
				histo.histo(), TOLERANCE );
		Assert.assertArrayEquals( new double[] {0, 0, 0, 3.0/5, 0, 0, 1.0/5}, 
				histo.histo(3), TOLERANCE );
		Assert.assertArrayEquals( new double[] {0, 1.0/5, 0, 0, 0, 0, 3.0/5, 0, 0, 1.0/5, 0, 0, 0}, 
				histo.histo(6), TOLERANCE );
	}
}
