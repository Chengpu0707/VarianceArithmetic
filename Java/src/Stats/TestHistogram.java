package Stats;

import org.junit.Assert;
import org.junit.Test;

public class TestHistogram {
	private static final double TOLERANCE = 3E-16;

	@Test
	public void testDefault() {
		final Histogram histo = new Histogram();
        Assert.assertEquals(5, histo.divids());
        Assert.assertEquals(3, histo.maxRange(), 1E-16);
        final double[] sExp = new double[31];
		Assert.assertTrue( null == histo.histo() );
		Assert.assertTrue( histo.accum( 0, 1 ) );
        sExp[15] = 1;
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( 0.01, 2 ) );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( -0.01, 3 ) );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		
		Assert.assertFalse( histo.accum( 5, 4 ) );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );

		Assert.assertTrue( histo.accum( -1, 5 ) );
		Assert.assertEquals( 0, histo.lower() );
        Assert.assertEquals( 1, histo.upper() );
        sExp[15] = 3.0/4;
        sExp[10] = 1.0/4;
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
		Assert.assertTrue( histo.accum( 0.295, 6 ) );
        sExp[15] = 3.0/5;
        sExp[10] = 1.0/5;
        sExp[16] = 1.0/5;
		Assert.assertEquals( 0, histo.lower() );
        Assert.assertEquals( 1, histo.upper() );
		Assert.assertArrayEquals( sExp, histo.histo(), TOLERANCE );
	}
}
