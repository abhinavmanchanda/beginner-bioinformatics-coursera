from functools import reduce

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    return pattern_count_with_mismatch(Text, Pattern, d)

def pattern_count_with_mismatch(text, pattern, mismatch_threshold):
    return sum(1 for i in range(len(text) - len(pattern) + 1) if HammingDistance(text[i:i + len(pattern)], pattern) <= mismatch_threshold)

def ApproximatePatternMatching(Text, Pattern, d):
    return list(i for i in range(len(Text) - len(Pattern) + 1) if HammingDistance(Text[i:i + len(Pattern)], Pattern) <= d)

def HammingDistance(p, q):
    return reduce(lambda a, b: a + b, list(map(lambda x, y: int(x != y), list(p), list(q))), 0)
