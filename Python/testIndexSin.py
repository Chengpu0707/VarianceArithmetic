import math
import unittest

from indexSin import IndexSin

class TestSin (unittest.TestCase):

    def test_neg_rem(self):
        '''
        assume order==3, so size==8
        '''
        # whe index is near 0
        self.assertEqual(0, 1 // 4)     #0 % 1
        self.assertEqual(1, 1 % 4)
        self.assertEqual(-1, -1 // 4)   #0 % -1
        self.assertEqual(3, -1 % 4)

        # whe index is near +pi/2
        self.assertEqual(0, 3 // 4)     #0 % 3
        self.assertEqual(3, 3 % 4)
        self.assertEqual(1, 5 // 4)     #0 % 3
        self.assertEqual(1, 5 % 4)

        # whe index is near +pi
        self.assertEqual(1, 7 // 4)     #0 % 1
        self.assertEqual(3, 7 % 4)
        self.assertEqual(2, 9 // 4)     #0 % -1
        self.assertEqual(1, 9 % 4)

        # whe index is near -pi/2
        self.assertEqual(-1, -3 // 4)   #0 % -3
        self.assertEqual(1, -3 % 4)
        self.assertEqual(-2, -5 // 4)   #0 % -3
        self.assertEqual(3, -5 % 4)

        # whe index is near -pi
        self.assertEqual(-2, -7 // 4)   #0 % -1
        self.assertEqual(1, -7 % 4)
        self.assertEqual(-3, -9 // 4)   #0 % 1
        self.assertEqual(3, -9 % 4)

    def test_get_sin_index(self):
        indexSin = IndexSin(3)
        self.assertEqual(8, indexSin.size())

        self.assertEqual( 0, indexSin.get_index(0))
        self.assertEqual( 1, indexSin.get_index(1))
        self.assertEqual( 2, indexSin.get_index(2))
        self.assertEqual( 3, indexSin.get_index(3))
        self.assertEqual( 4, indexSin.get_index(4))
        self.assertEqual( 3, indexSin.get_index(5))
        self.assertEqual( 2, indexSin.get_index(6))
        self.assertEqual( 1, indexSin.get_index(7))
        self.assertEqual( 0, indexSin.get_index(8))
        self.assertEqual(-1, indexSin.get_index(9))
        self.assertEqual(-2, indexSin.get_index(10))
        self.assertEqual(-3, indexSin.get_index(11))
        self.assertEqual(-4, indexSin.get_index(12))
        self.assertEqual(-3, indexSin.get_index(13))
        self.assertEqual(-2, indexSin.get_index(14))
        self.assertEqual(-1, indexSin.get_index(15))
        self.assertEqual( 0, indexSin.get_index(16))
        self.assertEqual( 1, indexSin.get_index(17))
        self.assertEqual( 2, indexSin.get_index(18))

        self.assertEqual(-1, indexSin.get_index(-1))
        self.assertEqual(-2, indexSin.get_index(-2))
        self.assertEqual(-3, indexSin.get_index(-3))
        self.assertEqual(-4, indexSin.get_index(-4))
        self.assertEqual(-3, indexSin.get_index(-5))
        self.assertEqual(-2, indexSin.get_index(-6))
        self.assertEqual(-1, indexSin.get_index(-7))
        self.assertEqual( 0, indexSin.get_index(-8))
        self.assertEqual( 1, indexSin.get_index(-9))
        self.assertEqual( 2, indexSin.get_index(-10))
        self.assertEqual( 3, indexSin.get_index(-11))
        self.assertEqual( 4, indexSin.get_index(-12))
        self.assertEqual( 3, indexSin.get_index(-13))
        self.assertEqual( 2, indexSin.get_index(-14))
        self.assertEqual( 1, indexSin.get_index(-15))
        self.assertEqual( 0, indexSin.get_index(-16))
        self.assertEqual(-1, indexSin.get_index(-17))
        self.assertEqual(-2, indexSin.get_index(-18))

    def test_sin(self):
        indexSin = IndexSin(3)
        self.assertEqual(8, indexSin.size())

        self.assertEqual(indexSin.sin(0), 0)
        self.assertEqual(indexSin.sin(1), math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(2), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(3), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(4), 1)
        self.assertEqual(indexSin.sin(5), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(6), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(7), math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(8), 0)

        self.assertEqual(indexSin.sin(-1), -math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(-2), -math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(-3), -math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(-4), -1)
        self.assertEqual(indexSin.sin(-5), -math.cos(1/8*math.pi))
        self.assertEqual(indexSin.sin(-6), -math.cos(2/8*math.pi))
        self.assertEqual(indexSin.sin(-7), -math.sin(1/8*math.pi))
        self.assertEqual(indexSin.sin(-8), 0)

    def test_sin_error(self):
        self.assertEqual(-math.sin(1/8*math.pi) - 9.992007221626409e-16, 
                         math.sin(1001/8*math.pi))
        self.assertEqual(math.ulp(math.sin(1/8*math.pi)), 5.551115123125783e-17)

    def test_cos(self):
        indexSin = IndexSin(3)
        self.assertEqual(8, indexSin.size())

        self.assertEqual(indexSin.cos(0), 1)
        self.assertEqual(indexSin.cos(1), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.cos(2), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.cos(3), math.sin(1/8*math.pi))
        self.assertEqual(indexSin.cos(4), 0)
        self.assertEqual(indexSin.cos(5), -math.sin(1/8*math.pi))
        self.assertEqual(indexSin.cos(6), -math.cos(2/8*math.pi))
        self.assertEqual(indexSin.cos(7), -math.cos(1/8*math.pi))
        self.assertEqual(indexSin.cos(8), -1)

        self.assertEqual(indexSin.cos(-1), math.cos(1/8*math.pi))
        self.assertEqual(indexSin.cos(-2), math.cos(2/8*math.pi))
        self.assertEqual(indexSin.cos(-3), math.sin(1/8*math.pi))

    def test_arc_sin(self):
        indexSin = IndexSin(3)
        self.assertEqual( 0, indexSin.arc_sin( 0))
        self.assertEqual( 1, indexSin.arc_sin( math.sin(math.pi*1/8)))
        self.assertEqual( 2, indexSin.arc_sin( math.sin(math.pi*2/8)))
        self.assertEqual( 3, indexSin.arc_sin( math.sin(math.pi*3/8)))
        self.assertEqual( 4, indexSin.arc_sin( 1))

        self.assertEqual(-1, indexSin.arc_sin(-math.sin(math.pi*1/8)))
        self.assertEqual(-2, indexSin.arc_sin(-math.sin(math.pi*2/8)))
        self.assertEqual(-3, indexSin.arc_sin(-math.sin(math.pi*3/8)))
        self.assertEqual(-4, indexSin.arc_sin(-1))

        self.assertEqual( 2.239966192595619, indexSin.arc_sin( math.sin(math.pi*1/3)))
        self.assertEqual(-2.239966192595619, indexSin.arc_sin(-math.sin(math.pi*1/3)))


if __name__ == '__main__':
    unittest.main()