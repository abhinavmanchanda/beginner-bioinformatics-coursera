from functools import reduce
from itertools import accumulate


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
    return reverse(dna_complement(Pattern))

def dna_complement(pattern):
    return ''.join(map(character_dna_complement, pattern))

def reverse(s):
    return s if len(s) == 0 else reverse(s[1:]) + s[0]

def character_dna_complement(character):
    return {"A": "T", "T": "A", "G": "C", "C": "G"}[character]

def half_string_symbol_count(text, symbol, index):
    length = len(text)
    extended_text = text + text[0:length // 2]
    return PatternCount(extended_text[index:index + (length // 2)], symbol)

def SymbolArray(genome, symbol):
    return {i: half_string_symbol_count(genome, symbol, i) for i in range(0, len(genome))}


def FasterSymbolArray(genome, symbol):
    initial_value = half_string_symbol_count(genome, symbol, 0)
    genome_length = len(genome)
    extended_genome = genome + genome[0:genome_length // 2]
    indices = list(range(1, genome_length))
    constructed_list = list(accumulate(indices,
               lambda result, index: update_half_array_symbol_count(extended_genome, index,
                                                                    genome_length, result, symbol), initial=initial_value))
    return dict(enumerate(constructed_list))

def update_half_array_symbol_count(extended_genome, index, original_genome_length, previous_value,
                                   symbol):
    return previous_value \
           - int(extended_genome[index - 1] == symbol) \
           + int(extended_genome[index + original_genome_length // 2 - 1] == symbol)
