def reverse_complement(pattern):
    complement = dna_complement(pattern)
    return reverse(complement)

def dna_complement(pattern):
    return ''.join(map(character_dna_complement, pattern))


def reverse(s):
    if len(s) == 0:
        return s
    else:
        return reverse(s[1:]) + s[0]


def character_dna_complement(pattern):
    dna_complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    return dna_complement[pattern]
