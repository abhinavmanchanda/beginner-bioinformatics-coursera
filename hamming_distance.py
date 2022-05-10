from functools import reduce

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    return len(ApproximatePatternMatching(Text, Pattern, d))

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(0, len(Text)+1-len(Pattern)):
        hamming_distance = HammingDistance(Text[i:i+len(Pattern)], Pattern)
        if is_approximate_match(hamming_distance, d):
            positions = positions + [i]
    return positions

def is_approximate_match(hamming_distance, threshold):
    return hamming_distance <= threshold


def HammingDistance(p, q):
    return reduce(lambda a, b: a + b, list(map(lambda x, y: int(x != y), list(p), list(q))), 0)
