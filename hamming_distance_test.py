import unittest
from hamming_distance import HammingDistance

class HammingDistanceTest(unittest.TestCase):
    def test_should_find_hamming_distance_between_two_strings(self):
        self.assertEqual(3, HammingDistance('GGGCCGTTGGT', 'GGACCGTTGAC'))
