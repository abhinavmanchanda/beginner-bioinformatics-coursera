from unittest import TestCase
from replication import count_occurence
from replication import frequency_map_for_kmers
from replication import most_frequent_substrings
from replication import ReverseComplement
from replication import half_string_symbol_count
from replication import SymbolArray
from replication import FasterSymbolArray
from replication import SkewArray
from replication import MinimumSkew

class Test(TestCase):
    def test_count_occurrence_of_a_pattern_in_a_text(self):
        self.assertEqual(2, count_occurence("GCGCG", "GCG"))

    def test_create_frequency_map_of_substrings_of_a_particular_length(self):
        expected_output = {'CGA': 1, 'GAT': 1, 'ATA': 3, 'TAT': 2, 'ATC': 1, 'TCC': 1, 'CCA': 1, 'CAT': 1, 'TAG': 1}
        self.assertEqual(expected_output, frequency_map_for_kmers("CGATATATCCATAG", 3))

    def test_find_most_frequent_substrings_of_a_given_length_in_a_string(self):
        expected_output = ['CATG', 'GCAT']
        self.assertCountEqual(expected_output, most_frequent_substrings('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4))

    def test_complement_and_reverse_to_find_pairing_of_DNA_strand(self):
        expected_output = 'ATACGTC'
        self.assertCountEqual(expected_output, ReverseComplement('GACGTAT'))

    def test_find_count_of_symbol_in_half_string_starting_at_index_in_circular_string(self):
        self.assertEqual(3, half_string_symbol_count('CACGCC', 'C', 4))

    def test_find_half_string_symbol_count_for_all_indices(self):
        expected_output = {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}
        self.assertEqual(expected_output, SymbolArray('AAAAGGGG', 'A'))

    def test_find_half_string_symbol_count_for_all_indices_in_o_n(self):
        expected_output = {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}
        self.assertEqual(expected_output, FasterSymbolArray('AAAAGGGG', 'A'))

    def test_create_skew_array_for_a_given_genome(self):
        expected_output = [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
        self.assertEqual(expected_output, SkewArray('CATGGGCATCGGCCATACGCC'))

    def test_find_min_skew_for_a_genome(self):
        expected_output = [11, 24]
        self.assertEqual(expected_output, MinimumSkew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'))