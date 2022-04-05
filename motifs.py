# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {}  # initializing the count dictionary
    motif_string_length = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for string_position in range(motif_string_length):
            count[symbol].append(0)

    number_of_motifs = len(Motifs)
    for current_motif in range(number_of_motifs):
        for string_position in range(motif_string_length):
            symbol = Motifs[current_motif][string_position]
            count[symbol][string_position] += 1
    return count


def Profile(Motifs):
    count = Count(Motifs)
    number_of_motifs = len(Motifs)
    profile = {}
    for symbol in "ACGT":
        profile[symbol] = motif_profile_for_symbol(count[symbol], number_of_motifs)
    return profile


def motif_profile_for_symbol(count_array_for_symbol, number_of_motifs):
    return list(map(lambda x: x / number_of_motifs, count_array_for_symbol))


def Consensus(Motifs):
    profile = Profile(Motifs)
    consensus_length = len(Motifs[0])
    symbols = "ACGT"
    consensus_profile = profile[symbols[0]]
    consensus = symbols[0] * consensus_length
    for symbol in symbols[1:len(symbols)]:
        symbol_profile = profile[symbol]
        for iterator in range(len(symbol_profile)):
            if symbol_profile[iterator] > consensus_profile[iterator]:
                consensus_profile[iterator] = symbol_profile[iterator]
                consensus = consensus[:iterator] + symbol + consensus[iterator + 1:]

    return consensus


from functools import reduce


def HammingDistance(p, q):
    plist = list(p)
    qlist = list(q)
    temp = list(map(lambda x, y: int(x != y), plist, qlist))
    return reduce(lambda a, b: a + b, temp, 0)


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