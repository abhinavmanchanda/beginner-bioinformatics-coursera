import unittest

from motifs import *


class MotifsTest(unittest.TestCase):

    def test_create_motif_count_dictionary_matrix(self):
        motifs = [
            "AACGTA",
            "CCCGTT",
            "CACCTT",
            "GGATTA",
            "TTCCGG"
        ]
        expected_output = {
            'A': [1, 2, 1, 0, 0, 2],
            'C': [2, 1, 4, 2, 0, 0],
            'G': [1, 1, 0, 2, 1, 1],
            'T': [1, 1, 0, 1, 4, 2]
        }

        self.assertEqual(expected_output, Count(motifs))

    def test_create_profile_matrix_for_motif_list(self):
        motifs = [
            "AACGTA",
            "CCCGTT",
            "CACCTT",
            "GGATTA",
            "TTCCGG"
        ]
        expected_output = {
            'A': [0.2, 0.4, 0.2, 0, 0, 0.4],
            'C': [0.4, 0.2, 0.8, 0.4, 0, 0],
            'G': [0.2, 0.2, 0, 0.4, 0.2, 0.2],
            'T': [0.2, 0.2, 0, 0.2, 0.8, 0.4]
        }

        self.assertEqual(expected_output, Profile(motifs))

    def test_find_consensus_symbol_for_a_particular_index_in_list_of_motifs(self):
        motifs = ["AACGTA",
                  "CCCGTT",
                  "CACCTT",
                  "GGATTA",
                  "TTCCGG"]
        self.assertEqual("C", consensus_symbol_at_index(Count(motifs), 2, "ACGT"))

    def test_find_consensus_string_for_a_given_list_of_motifs(self):
        motifs = ["AACGTA",
                  "CCCGTT",
                  "CACCTT",
                  "GGATTA",
                  "TTCCGG"]
        actual_output = Consensus(motifs)
        possible_expected_outputs = ["CACCTA", "CACGTT", "CACCTT", "CACGTA"]
        self.assertTrue(actual_output in possible_expected_outputs)

    def test_consensus_string_score_for_a_list_of_motifs(self):
        self.assertEqual(14, Score(["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]))

    def test_probability_of_motif_from_profile_matrix(self):
        expected = 0.22222
        actual = Pr("AGA", {'A': [0.33, 0, 1], 'T': [0.67, 0, 0], 'C': [0, 0.33, 0], 'G': [0, 0.67, 0]})
        self.assertAlmostEqual(expected, actual, 2)

    def test_find_most_probable_kmer_in_string_from_a_profile_matrix(self):
        expected = "TGA"
        actual = ProfileMostProbableKmer("AGATGACCA", 3,
                        {'A': [0.33, 0, 1], 'T': [0.67, 0, 0], 'C': [0, 0.33, 0], 'G': [0, 0.67, 0]})
        self.assertEqual(expected, actual)

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

    def test_find_best_motif_out_of_multiple_randomized_runs(self):
        dna = [
            'GGCGTTCAGGCA',
            'AAGAATCAGTCA',
            'CAAGGAGTTCGC',
            'CACGTCAATCAC',
            'CAATAATATTCG'
        ]
        self.assertEqual(Score(['CAG', 'CAG', 'CAA', 'CAA', 'CAA']), Score(best_randomised_motifs(dna, 3, 5)))

    def test_normalize_probabilities(self):
        self.assertEqual({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
                         Normalize({'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1})
                         )

    def test_create_probability_ranges_from_probabilities(self):
        self.assertEqual(
            {'A': {'lower': 0, 'upper': 0.25},
             'B': {'lower': 0.25, 'upper': 4.25},
             'C': {'lower': 4.25, 'upper': 16.25},
             'D': {'lower': 16.25, 'upper': 27.36}
             },
            key_vs_range({'A': 0.25, 'B': 4, 'C': 12, 'D': 11.11}))

    def test_weighted_die_throw(self):
        output_count = {'x': 0, 'y': 0, 'abc': 0}
        for i in range(10000):
            current_output = WeightedDie({'x': 0.33, 'y': 0.5, 'abc': .17})
            output_count[current_output] = output_count[current_output] + 1
        self.assertTrue(4400 < output_count['y'] < 5500)
        self.assertTrue(2800 < output_count['x'] < 3800)
        self.assertTrue(1200 < output_count['abc'] < 2200)

    def test_generate_weighted_probability_motif_based_on_profile(self):
        text = 'AAACCCAAACCC'
        profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
        output_count = dict()
        for i in range(10000):
            current_output = ProfileGeneratedString(text, profile, 2)
            output_count[current_output] = output_count.get(current_output, 0) + 1
        self.assertTrue(1500 < output_count['AA'] < 2500)
        self.assertTrue(3600 < output_count['AC'] < 4600)
        self.assertTrue(2000 < output_count['CC'] < 3000)
        self.assertTrue(800 < output_count['CA'] < 1800)

    def test_find_best_motifs_using_gibbs_sampling(self):
        Dna = [
            "GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
            "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
            "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
            "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
            "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
            "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
            "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
            "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
            "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
            "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"
        ]
        kmer_size = 8
        number_of_dna_strings = 5
        number_of_iterations = 100

        best_score = 100
        for i in range(100):
            actual = GibbsSampler(Dna, kmer_size, number_of_dna_strings, number_of_iterations)
            if Score(actual) < best_score :
                best_score = Score(actual)

        self.assertTrue(best_score < 40)

    def test_generate_next_gibbs_sampler_motif_for_index(self):
        dna = [
            'TTACCTTAAC',
            'GATGTCTGTC',
            'CCGGCGTTAG',
            'CACTAACGAG',
            'CGTCAGAGGT'
        ]
        current_motifs = [
            'ACCT',
            'GTCT',
            'GCGT',
            'ACTA',
            'AGGT'
        ]

        acga_count = 0
        for i in range(1000) :
            output = generate_new_motifs(dna_strings=dna, kmer_length=4, current_motifs=current_motifs, index=3)
            if output[3] == 'ACGA':
                acga_count += 1
        self.assertTrue(400 < acga_count < 430)

if __name__ == '__main__':
    unittest.main()
