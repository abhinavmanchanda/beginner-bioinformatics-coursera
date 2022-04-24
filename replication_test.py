from unittest import TestCase
from replication import count_occurence
from replication import frequency_map_for_kmers
from replication import most_frequent_substrings
from replication import ReverseComplement

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
