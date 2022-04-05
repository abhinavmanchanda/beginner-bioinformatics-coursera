from functools import reduce


def SkewArray(Genome):
    nucleotide_to_value = {"A":0, "T":0, "G":1, "C":-1}
    output = [0]
    for i in range(1, len(Genome)+1):
        output += [output[i-1] + nucleotide_to_value[Genome[i-1]]]
    return output

print(SkewArray("CATGGGCATCGGCCATACGCC"))