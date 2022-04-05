# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    length = len(Genome)
    extended_genome = Genome + Genome[0:length//2]
    array[0] = extended_genome[0:length//2].count(symbol)
    for i in range(1, length):
        array[i] = array[i-1] - extended_genome[i-1:i].count(symbol) + extended_genome[i+length//2-1:i+length//2].count(symbol)

    return array

file = open('./e_coli.txt','r')
ecoli_genome = file.read()
file.close()
array = FasterSymbolArray(ecoli_genome, 'C')

import matplotlib.pyplot as plt
plt.plot(*zip(*sorted(array.items())))
plt.show()