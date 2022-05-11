import unittest
from hamming_distance import HammingDistance
from hamming_distance import ApproximatePatternCount
from hamming_distance import ApproximatePatternMatching

class HammingDistanceTest(unittest.TestCase):
    def test_find_hamming_distance_between_two_strings(self):
        self.assertEqual(3, HammingDistance('GGGCCGTTGGT', 'GGACCGTTGAC'))

    def test_find_approximate_pattern_count_with_mismatch_threshold(self):
        self.assertEqual(4, ApproximatePatternCount('GAGG', 'TTTAGAGCCTTCAGAGG', 2))

    def test_find_indices_of_approximate_pattern_matches_with_given_threshold(self):
        pattern = 'ATTCTGGA'
        dna = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        mismatch_threshold = 3
        self.assertEqual([6, 7, 26, 27], ApproximatePatternMatching(dna, pattern, mismatch_threshold))