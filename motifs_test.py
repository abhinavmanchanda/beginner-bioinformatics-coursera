import unittest
from motifs import GreedyMotifSearch
from motifs import CountWithPseudocounts
from motifs import ProfileWithPseudocounts
from motifs import GreedyMotifSearchWithPseudocounts
from motifs import Motifs


class MotifsTest(unittest.TestCase):
    def test_should_greedily_search_for_motif(self):
        k = 3
        t = 5
        Dna = [
            "GGCGTTCAGGCA",
            "AAGAATCAGTCA",
            "CAAGGAGTTCGC",
            "CACGTCAATCAC",
            "CAATAATATTCG"
        ]
        expected_output = [
            "CAG",
            "CAG",
            "CAA",
            "CAA",
            "CAA"
        ]
        self.assertEqual(expected_output, GreedyMotifSearch(Dna, k, t))

    def test_add_pseudocount_of_one_to_each_motif_count(self):
        motifs = [
            "AACGTA",
            "CCCGTT",
            "CACCTT",
            "GGATTA",
            "TTCCGG"
        ]
        expected_output = {'A': [2, 3, 2, 1, 1, 3], 'C': [3, 2, 5, 3, 1, 1], 'G': [2, 2, 1, 3, 2, 2],
                           'T': [2, 2, 1, 2, 5, 3]}

        self.assertEqual(expected_output, CountWithPseudocounts(motifs))

    def test_create_profile_with_pseudocounts(self):
        motifs = [
            "AACGTA",
            "CCCGTT",
            "CACCTT",
            "GGATTA",
            "TTCCGG"
        ]
        expected_output = {
            'A': [0.2222222222222222, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111,
                  0.3333333333333333],
            'C': [0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.1111111111111111,
                  0.1111111111111111],
            'G': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222,
                  0.2222222222222222],
            'T': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556,
                  0.3333333333333333]}
        self.assertEqual(expected_output, ProfileWithPseudocounts(motifs))

    def test_greedy_motif_search_with_pseudocount(self):
        k = 3
        t = 5
        Dna = [
            "GGCGTTCAGGCA",
            "AAGAATCAGTCA",
            "CAAGGAGTTCGC",
            "CACGTCAATCAC",
            "CAATAATATTCG"
        ]
        expected_output = [
            "TTC",
            "ATC",
            "TTC",
            "ATC",
            "TTC"
        ]
        self.assertEqual(expected_output, GreedyMotifSearchWithPseudocounts(Dna, k, t))

    def test_find_most_probable_motifs_from_dna_given_profile_matrix(self):
        profile = {'A': [0.8, 0.0, 0.0, 0.2],
                   'C': [0.0, 0.6, 0.2, 0.0],
                   'G': [0.2, 0.2, 0.8, 0.0],
                   'T': [0.0, 0.2, 0.0, 0.8]}
        dna = ["TTACCTTAAC",
               "GATGTCTGTC",
               "ACGGCGTTAG",
               "CCCTAACGAG",
               "CGTCAGAGGT"
               ]

        expected_output = [
            "ACCT",
            "ATGT",
            "GCGT",
            "ACGA",
            "AGGT"
        ]
        self.assertEqual(expected_output, Motifs(profile, dna))


if __name__ == '__main__':
    unittest.main()
