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

from functools import reduce

def HammingDistance(p, q):
    return sum(list(map(lambda x, y: int(x != y), list(p), list(q))))

def Score(Motifs):
    consensus = Consensus(Motifs)
    diff_list = map(lambda motif: HammingDistance(motif, consensus), Motifs)
    return reduce(lambda a, b: a + b, diff_list)


def Pr(Text, Profile):
    probability = 1
    for i in range(len(Text)):
        probability *= Profile[Text[i]][i]
    return probability


def ProfileMostProbableKmer(text, k, profile):
    probabilities = [0] * (len(text) - k + 1)
    for iterator in range(len(text) - k + 1):
        kmer = text[iterator:iterator + k]
        probabilities[iterator] = Pr(kmer, profile)
    probable_index = probabilities.index(max(probabilities))
    return text[probable_index:probable_index + k]


def GreedyMotifSearch(Dna, k, t):
    best_motifs = list(map(lambda current_dna: current_dna[0:k], Dna))
    dna_string_length = len(Dna[0])
    for i in range(dna_string_length - k + 1):
        current_motifs = [""] * t
        current_motifs[0] = Dna[0][i:i + k]
        for j in range(1, t):
            profile = Profile(current_motifs[0:j])
            current_motifs[j] = ProfileMostProbableKmer(Dna[j], k, profile)
        if Score(current_motifs) < Score(best_motifs):
            best_motifs = current_motifs
    return best_motifs


def CountWithPseudocounts(Motifs):
    motifs_count = Count(Motifs)
    return {key: add_pseudocount_toarray(value) for key, value in motifs_count.items()}


def add_pseudocount_toarray(motif_count_array):
    return list(map(lambda x: x + 1, motif_count_array))


# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    motifs_pseudocounts = CountWithPseudocounts(Motifs)
    divisor = len(Motifs) + 4
    return {key: list(map(lambda x: x / divisor, value)) for key, value in motifs_pseudocounts.items()}


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    best_motifs = list(map(lambda current_dna: current_dna[0:k], Dna))
    dna_string_length = len(Dna[0])
    for i in range(dna_string_length - k + 1):
        current_motifs = [""] * t
        current_motifs[0] = Dna[0][i:i + k]
        for j in range(1, t):
            profile = ProfileWithPseudocounts(current_motifs[0:j])
            current_motifs[j] = ProfileMostProbableKmer(Dna[j], k, profile)
        if Score(current_motifs) < Score(best_motifs):
            best_motifs = current_motifs
    return best_motifs


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
