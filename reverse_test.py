import unittest
from reverse import reverse_complement


class ReverseTest(unittest.TestCase):
    def test_reverse_and_complement_empty_string(self):
        self.assertEqual("", reverse_complement(""))

    def test_reverse_complement_single_letter_A(self):
        self.assertEqual("T", reverse_complement("A"))

    def test_reverse_complement_single_letters_CGT(self):
        self.assertEqual("C", reverse_complement("G"))
        self.assertEqual("G", reverse_complement("C"))
        self.assertEqual("A", reverse_complement("T"))

    def test_reverse_complement_two_characters(self):
        self.assertEqual("A", reverse_complement("GT"))


if __name__ == '__main__':
    unittest.main()
