package Stats;

import org.junit.Assert;
import org.junit.Test;

public class TestStat extends Stat {
	
	private static final double TOLERANCE = 1E-15;

	@Test
	public void test() {
		Assert.assertEquals( 0, count() );
		Assert.assertTrue( Double.isNaN( avg() ) );
		Assert.assertEquals( Double.MAX_VALUE, min(), TOLERANCE );
		Assert.assertEquals( -Double.MAX_VALUE, max(), TOLERANCE );

		Assert.assertFalse( accum( Double.NaN, 1 ) );
		Assert.assertFalse( accum( Double.POSITIVE_INFINITY, 2 ) );
		Assert.assertFalse( accum( Double.NEGATIVE_INFINITY, 3 ) );
		
		Assert.assertEquals( 0, count() );
		Assert.assertTrue( Double.isNaN( avg() ) );
		Assert.assertEquals( Double.MAX_VALUE, min(), TOLERANCE );
		Assert.assertEquals( -Double.MAX_VALUE, max(), TOLERANCE );
		
		Assert.assertTrue( accum( 0, 1 ) );
		Assert.assertEquals( 1, count() );
		Assert.assertEquals( 0, avg(), TOLERANCE );
		Assert.assertEquals( 0, var(), TOLERANCE );
		Assert.assertEquals( 0, dev(), TOLERANCE );
		Assert.assertEquals( 0, min(), TOLERANCE );
		Assert.assertEquals( 0, max(), TOLERANCE );

		Assert.assertTrue( accum( 1, 2 ) );
		Assert.assertTrue( accum( -1, 3 ) );
		Assert.assertEquals( 3, count() );
		Assert.assertEquals( 0, avg(), TOLERANCE );
		Assert.assertEquals( 2.0/3, var(), TOLERANCE );
		Assert.assertEquals( Math.sqrt(2.0/3), dev(), TOLERANCE );
		Assert.assertEquals( -1, min(), TOLERANCE );
		Assert.assertEquals( 1, max(), TOLERANCE );
	}

}
