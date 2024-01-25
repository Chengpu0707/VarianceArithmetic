import math
import unittest

from indexSin import IndexSin

class TestSin (unittest.TestCase):

    def test_neg_rem(self):
        '''
        assume order==3, so size==8
        '''
        # whe index is near 0
        self.assertEqual(0, 1 // 4)
        self.assertEqual(1, 1 % 4)

        self.assertEqual(-1, -1 // 4)
        self.assertEqual(3, -1 % 4)

        # whe index is near +pi/2
        self.assertEqual(0, 3 // 4)
        self.assertEqual(3, 3 % 4)
        self.assertEqual(1, 5 // 4)
        self.assertEqual(1, 5 % 4)

        # whe index is near -pi/2
        self.assertEqual(-1, -3 // 4)
        self.assertEqual(1, -3 % 4)
        self.assertEqual(-2, -5 // 4)
        self.assertEqual(3, -5 % 4)

        # whe index is near +pi
        self.assertEqual(1, 7 // 4)
        self.assertEqual(3, 7 % 4)
        self.assertEqual(2, 9 // 4)
        self.assertEqual(1, 9 % 4)

        # whe index is near -pi
        self.assertEqual(-2, -7 // 4)
        self.assertEqual(1, -7 % 4)
        self.assertEqual(-3, -9 // 4)
        self.assertEqual(3, -9 % 4)

    def test_get_sin_index(self):
        index_sin = IndexSin(3)
        self.assertEqual(8, index_sin.size())

        self.assertEqual( 0, index_sin._get_sin_index(0))
        self.assertEqual( 1, index_sin._get_sin_index(1))
        self.assertEqual( 2, index_sin._get_sin_index(2))
        self.assertEqual( 3, index_sin._get_sin_index(3))
        self.assertEqual( 4, index_sin._get_sin_index(4))
        self.assertEqual( 3, index_sin._get_sin_index(5))
        self.assertEqual( 2, index_sin._get_sin_index(6))
        self.assertEqual( 1, index_sin._get_sin_index(7))
        self.assertEqual( 0, index_sin._get_sin_index(8))
        self.assertEqual(-1, index_sin._get_sin_index(9))
        self.assertEqual(-2, index_sin._get_sin_index(10))
        self.assertEqual(-3, index_sin._get_sin_index(11))
        self.assertEqual(-4, index_sin._get_sin_index(12))
        self.assertEqual(-3, index_sin._get_sin_index(13))
        self.assertEqual(-2, index_sin._get_sin_index(14))
        self.assertEqual(-1, index_sin._get_sin_index(15))
        self.assertEqual( 0, index_sin._get_sin_index(16))
        self.assertEqual( 1, index_sin._get_sin_index(17))
        self.assertEqual( 2, index_sin._get_sin_index(18))

        self.assertEqual(-1, index_sin._get_sin_index(-1))
        self.assertEqual(-2, index_sin._get_sin_index(-2))
        self.assertEqual(-3, index_sin._get_sin_index(-3))
        self.assertEqual(-4, index_sin._get_sin_index(-4))
        self.assertEqual(-3, index_sin._get_sin_index(-5))
        self.assertEqual(-2, index_sin._get_sin_index(-6))
        self.assertEqual(-1, index_sin._get_sin_index(-7))
        self.assertEqual( 0, index_sin._get_sin_index(-8))
        self.assertEqual( 1, index_sin._get_sin_index(-9))
        self.assertEqual( 2, index_sin._get_sin_index(-10))
        self.assertEqual( 3, index_sin._get_sin_index(-11))
        self.assertEqual( 4, index_sin._get_sin_index(-12))
        self.assertEqual( 3, index_sin._get_sin_index(-13))
        self.assertEqual( 2, index_sin._get_sin_index(-14))
        self.assertEqual( 1, index_sin._get_sin_index(-15))
        self.assertEqual( 0, index_sin._get_sin_index(-16))
        self.assertEqual(-1, index_sin._get_sin_index(-17))
        self.assertEqual(-2, index_sin._get_sin_index(-18))

    def test_get_cos_index(self):
        index_sin = IndexSin(3)
        self.assertEqual(8, index_sin.size())

        self.assertEqual( 4, index_sin._get_cos_index(0))
        self.assertEqual( 3, index_sin._get_cos_index(1))
        self.assertEqual( 2, index_sin._get_cos_index(2))
        self.assertEqual( 1, index_sin._get_cos_index(3))
        self.assertEqual( 0, index_sin._get_cos_index(4))
        self.assertEqual(-1, index_sin._get_cos_index(5))
        self.assertEqual(-2, index_sin._get_cos_index(6))
        self.assertEqual(-3, index_sin._get_cos_index(7))
        self.assertEqual(-4, index_sin._get_cos_index(8))
        self.assertEqual(-3, index_sin._get_cos_index(9))
        self.assertEqual(-2, index_sin._get_cos_index(10))
        self.assertEqual(-1, index_sin._get_cos_index(11))
        self.assertEqual( 0, index_sin._get_cos_index(12))
        self.assertEqual( 1, index_sin._get_cos_index(13))
        self.assertEqual( 2, index_sin._get_cos_index(14))
        self.assertEqual( 3, index_sin._get_cos_index(15))
        self.assertEqual( 4, index_sin._get_cos_index(16))
        self.assertEqual( 3, index_sin._get_cos_index(17))
        self.assertEqual( 2, index_sin._get_cos_index(18))

        self.assertEqual( 3, index_sin._get_cos_index(-1))
        self.assertEqual( 2, index_sin._get_cos_index(-2))
        self.assertEqual( 1, index_sin._get_cos_index(-3))
        self.assertEqual( 0, index_sin._get_cos_index(-4))
        self.assertEqual(-1, index_sin._get_cos_index(-5))
        self.assertEqual(-2, index_sin._get_cos_index(-6))
        self.assertEqual(-3, index_sin._get_cos_index(-7))
        self.assertEqual(-4, index_sin._get_cos_index(-8))
        self.assertEqual(-3, index_sin._get_cos_index(-9))
        self.assertEqual(-2, index_sin._get_cos_index(-10))
        self.assertEqual(-1, index_sin._get_cos_index(-11))
        self.assertEqual( 0, index_sin._get_cos_index(-12))
        self.assertEqual( 1, index_sin._get_cos_index(-13))
        self.assertEqual( 2, index_sin._get_cos_index(-14))
        self.assertEqual( 3, index_sin._get_cos_index(-15))
        self.assertEqual( 4, index_sin._get_cos_index(-16))
        self.assertEqual( 3, index_sin._get_cos_index(-17))
        self.assertEqual( 2, index_sin._get_cos_index(-18))

    def test_sin(self):
        index_sin = IndexSin(3)
        self.assertEqual(8, index_sin.size())

        self.assertEqual(0, index_sin.sin(0))
        self.assertEqual(math.sin(1/8*math.pi), index_sin.sin(1))
        self.assertEqual(math.cos(2/8*math.pi), index_sin.sin(2))
        self.assertEqual(math.cos(1/8*math.pi), index_sin.sin(3))
        self.assertEqual(1, index_sin.sin(4))
        self.assertEqual(math.cos(1/8*math.pi), index_sin.sin(5))
        self.assertEqual(math.cos(2/8*math.pi), index_sin.sin(6))
        self.assertEqual(math.sin(1/8*math.pi), index_sin.sin(7))
        self.assertEqual(0, index_sin.sin(8))

        self.assertEqual(-math.sin(1/8*math.pi), index_sin.sin(-1))
        self.assertEqual(-math.cos(2/8*math.pi), index_sin.sin(-2))
        self.assertEqual(-math.cos(1/8*math.pi), index_sin.sin(-3))
        self.assertEqual(-1, index_sin.sin(-4))
        self.assertEqual(-math.cos(1/8*math.pi), index_sin.sin(-5))
        self.assertEqual(-math.cos(2/8*math.pi), index_sin.sin(-6))
        self.assertEqual(-math.sin(1/8*math.pi), index_sin.sin(-7))
        self.assertEqual(0, index_sin.sin(8))

        self.assertEqual(-math.sin(1/8*math.pi) - 9.992007221626409e-16, math.sin(1001/8*math.pi))
        self.assertEqual(-math.ulp(math.sin(1/8*math.pi)) + 5.551115123125783e-17, 0)

    def test_cos(self):
        index_sin = IndexSin(3)
        self.assertEqual(8, index_sin.size())

        self.assertEqual(1, index_sin.cos(0))
        self.assertEqual( math.cos(1/8*math.pi), index_sin.cos(1))
        self.assertEqual( math.cos(2/8*math.pi), index_sin.cos(2))
        self.assertEqual( math.sin(1/8*math.pi), index_sin.cos(3))
        self.assertEqual(0, index_sin.cos(4))
        self.assertEqual(-math.sin(1/8*math.pi), index_sin.cos(5))
        self.assertEqual(-math.cos(2/8*math.pi), index_sin.cos(6))
        self.assertEqual(-math.cos(1/8*math.pi), index_sin.cos(7))
        self.assertEqual(-1, index_sin.cos(8))

        self.assertEqual( math.cos(1/8*math.pi), index_sin.cos(-1))
        self.assertEqual( math.cos(2/8*math.pi), index_sin.cos(-2))
        self.assertEqual( math.sin(1/8*math.pi), index_sin.cos(-3))


if __name__ == '__main__':
    unittest.main()