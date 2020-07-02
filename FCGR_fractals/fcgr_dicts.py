import collections
from os import getcwd
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import math


loc = getcwd() + "/data/Test1/Anelloviridae/Anelloviridae_199.fasta"
f = open(loc)
s1 = f.read()
data = "".join(s1.split("\n")[1:])
data

#The N is used to remove a certain base pair A,T,C,G if needed
#finding the number of each unique k-mer in the sequence, for example
#it finds the number of time the kmer CAT appears in the sequence, and continues with each 3-mer
def count_kmers(sequence, k):
    d = collections.defaultdict(int)
    for i in range(len(data)-(k-1)):
        d[sequence[i:i+k]] +=1
    for key in d.keys():
        if "N" in key:
            del d[key]
    return d


# getting the count of a specific kmer, dividing by (length of sequence - length of kmer +1)
def probabilities(kmer_count, k):
        probabilities = collections.defaultdict(float)
        N = len(data)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - k + 1)
        return probabilities

def chaos_game_representation(probabilities, k):
        array_size = int(math.sqrt(4**k))
        chaos = []
        for i in range(array_size):
            chaos.append([0]*array_size)
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        for key, value in probabilities.items():
            for char in key:
                if char == "T":
                    posx += maxx / 2
                    posx = int(posx)
                elif char == "C":
                    posy += maxy / 2
                    posy = int(posy)
                elif char == "G":
                    posx += maxx / 2
                    posx = int(posx)
                    posy += maxy / 2
                    posy = int(posy)
                maxx = maxx / 2
                maxy /= 2
            chaos[posy-1][posx-1] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1

        return chaos

#each box represents a k-mer
f3 = count_kmers(data, 3)
len(data)
f3
f3_prob = probabilities(f3, 3)
chaos_k3 = chaos_game_representation(f3_prob, 3)
plt.title('Chaos game representation for 3-mers')
plt.imshow(chaos_k3, interpolation='nearest', cmap=cm.gray_r)
plt.show()

f4 = count_kmers(data, 4)
f4_prob = probabilities(f4, 4)
chaos_k4 = chaos_game_representation(f4_prob, 4)
plt.title('Chaos game representation for 4-mers')
plt.imshow(chaos_k4, interpolation='nearest', cmap=cm.gray_r)
plt.show()



f7 = count_kmers(data, 7)
f7_prob = probabilities(f7, 7)
chaos_k7 = chaos_game_representation(f7_prob, 7)
plt.title('Chaos game representation for 7-mers')
plt.imshow(chaos_k7, interpolation='nearest', cmap=cm.gray_r)
plt.show()
