from functools import reduce


def PatternCount(Text, Pattern):
    return count_occurence(Text, Pattern)


def count_occurence(text, pattern):
    return sum(1 for i in range(len(text) - len(pattern) + 1) if text[i:i + len(pattern)] == pattern)


def combine_dictionary_counts(a, b):
    return dict(list(a.items()) + list(b.items()) + [(k, a[k] + b[k]) for k in set(b) & set(a)])

def frequency_map_for_kmers(text, length):
    return reduce(combine_dictionary_counts, [{text[i:i + length]: 1} for i in range(len(text) - length + 1)])

def most_frequent_substrings(text, length):
    frequency_map = frequency_map_for_kmers(text, length)
    max_frequency = max(frequency for kmer, frequency in frequency_map.items())
    return [kmer for kmer, frequency in frequency_map.items() if frequency == max_frequency]

def FrequencyMap(Text, k):
    return frequency_map_for_kmers(Text, k)

def FrequentWords(Text, k):
    return most_frequent_substrings(Text, k)

def ReverseComplement(Pattern):
    complement = dna_complement(Pattern)
    return reverse(complement)

def dna_complement(pattern):
    return ''.join(map(character_dna_complement, pattern))


def reverse(s):
    return s if len(s) == 0 else reverse(s[1:]) + s[0]

def character_dna_complement(character):
    return {"A": "T", "T": "A", "G": "C", "C": "G"}[character]
