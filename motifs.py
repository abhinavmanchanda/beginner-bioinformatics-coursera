def Count(Motifs):
    return {symbol: count_array(symbol, Motifs) for symbol in "ACGT"}

def count_array(symbol, motifs):
   individual_count_arrays = [symbol_count_array_for_single_motif(symbol, motif) for motif in motifs]
   return [sum(elements) for elements in zip(*individual_count_arrays)]

def symbol_count_array_for_single_motif(symbol, motif):
    return [1 if current_symbol == symbol else 0 for current_symbol in motif]

def Profile(Motifs):
    count = Count(Motifs)
    number_of_motifs = len(Motifs)
    return {symbol: motif_profile_for_symbol(count[symbol], number_of_motifs) for symbol in "AGCT"}


def motif_profile_for_symbol(count_array_for_symbol, number_of_motifs):
    return list(map(lambda x: x / number_of_motifs, count_array_for_symbol))

def Consensus(Motifs):
    count_matrix = Count(Motifs)
    symbols = "ACGT"
    string_length = len(Motifs[0])
    consensus_list = [consensus_symbol_at_index(count_matrix, i, symbols) for i in range(string_length)]
    return ''.join(consensus_list)

def consensus_symbol_at_index(count_matrix, index, symbols):
    count_to_symbols_tuple = {count_matrix[symbol][index]: symbol for symbol in symbols}
    max_count = max(count_to_symbols_tuple.keys())
    return count_to_symbols_tuple[max_count]

def HammingDistance(p, q):
    return sum(list(map(lambda x, y: int(x != y), list(p), list(q))))

def Score(Motifs):
    return sum(map(lambda motif: HammingDistance(motif, Consensus(Motifs)), Motifs))

from functools import reduce

def probability_of_generation(motif, profile_matrix):
    return reduce(lambda x, y: x * y, [profile_matrix[motif[i]][i] for i in range(len(motif))], 1)

def Pr(Text, Profile):
    return probability_of_generation(Text, Profile)

def ProfileMostProbableKmer(text, k, profile):
    kmers = [text[iterator:iterator + k] for iterator in range(len(text) - k + 1)]
    probabilities = [probability_of_generation(kmer, profile) for kmer in kmers]
    return kmers[probabilities.index(max(probabilities))]

def GreedyMotifSearch(Dna, k, t):
    return greedy_search_with_custom_profile_function(Dna, k, t, Profile)

def greedy_search_with_custom_profile_function(Dna, k, t, profile_function):
    motif_combinations = [best_motifs_for_given_iteration(Dna, k, i, profile_function) for i in range(len(Dna[0]) - k + 1)]
    motif_scores = [Score(motifs) for motifs in motif_combinations]
    return motif_combinations[motif_scores.index(min(motif_scores))]

def best_motifs_for_given_iteration(dna, substring_length, index, profile_function):
    substring = dna[0][index: index + substring_length]
    profile_matrix = profile_function([substring])
    return recursive_compute_best_motifs(dna, substring_length, [substring], profile_matrix, 1, profile_function)

def recursive_compute_best_motifs(dna, substring_length, previous_motifs, profile_matrix, row_index, profile_function):
    if row_index == len(dna):
        return previous_motifs
    motif_for_row_index = ProfileMostProbableKmer(dna[row_index], substring_length, profile_matrix)
    current_motifs = previous_motifs + [motif_for_row_index]
    return recursive_compute_best_motifs(dna, substring_length, current_motifs,
                                         profile_function(current_motifs), row_index + 1, profile_function)


def CountWithPseudocounts(Motifs):
    motifs_count = Count(Motifs)
    return {key: add_pseudocount_toarray(value) for key, value in motifs_count.items()}

def add_pseudocount_toarray(motif_count_array):
    return list(map(lambda x: x + 1, motif_count_array))

def ProfileWithPseudocounts(Motifs):
    motifs_pseudocounts = CountWithPseudocounts(Motifs)
    divisor = len(Motifs) + 4
    return {key: list(map(lambda x: x / divisor, value)) for key, value in motifs_pseudocounts.items()}

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    return greedy_search_with_custom_profile_function(Dna, k, t, ProfileWithPseudocounts)

def Motifs(Profile, Dna):
    return list(map(lambda text: ProfileMostProbableKmer(text, len(Profile[next(iter(Profile))]), Profile), Dna))

import random

def RandomMotifs(Dna, k, t):
    return list(map(lambda text: select_random_motif(text, k), Dna))


def select_random_motif(dna_text, motif_length):
    position = random.randint(0, len(dna_text) - motif_length)
    return dna_text[position: position + motif_length]


def RandomizedMotifSearch(Dna, k, t):
    random_motifs = RandomMotifs(Dna, k, t)
    return converge_to_optimum_motifs(random_motifs, Dna)


def converge_to_optimum_motifs(current_motifs, dna):
    newly_computed_motifs = Motifs(ProfileWithPseudocounts(current_motifs), dna)
    if Score(newly_computed_motifs) == Score(current_motifs):
        return current_motifs
    return converge_to_optimum_motifs(newly_computed_motifs, dna)


def best_randomised_motifs(dna_array, motif_length, runs):
    best_of_current_run = RandomizedMotifSearch(dna_array, motif_length, len(dna_array))
    if runs == 1: return best_of_current_run
    best_of_other_runs = best_randomised_motifs(dna_array, motif_length, runs - 1)
    return best_of_current_run if Score(best_of_current_run) < Score(best_of_other_runs) else best_of_other_runs

def Normalize(Probabilities):
    sum_of_probabilities = sum(Probabilities.values())
    return {k: v/sum_of_probabilities for k, v in Probabilities.items()}

def WeightedDie(Probabilities):
    key_to_range = key_vs_range(Probabilities)
    random_fraction = random.uniform(0, 1)
    return next(key for key, value in key_to_range.items() if value['lower'] < random_fraction <= value['upper'])


def key_vs_range(Probabilities):
    key_list = Probabilities.keys()
    value_list = map(lambda key: Probabilities[key], key_list)
    upper_bound_list = reduce(lambda result, element: result + [result[-1] + element], value_list, [0])[1:]
    lower_bound_list = [0] + upper_bound_list[:-1]
    lower_and_upper_bounds = list(
        map(lambda lower, upper: {'lower': lower, 'upper': upper}, lower_bound_list, upper_bound_list))
    return dict(zip(key_list, lower_and_upper_bounds))

def ProfileGeneratedString(Text, profile, k):
    text_length = len(Text)
    probabilities = {Text[index:index+k]: Pr(Text[index:index+k], profile) for index in range(text_length - k + 1)}
    weighted_probabilities = Normalize(probabilities)
    return WeightedDie(weighted_probabilities)
