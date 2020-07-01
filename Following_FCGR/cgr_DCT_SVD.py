import collections
import numpy as np
import pylab
import math
from os import getcwd
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.fftpack import dct, idct
import pandas as pd
from PIL import Image
import cv2

#_______________________________________________________________________________
# Finding FCGR

# Extracting Sequence from File
loc = getcwd() + "/data/Test1/Anelloviridae/Anelloviridae_199.fasta"
f = open(loc)
s1 = f.read()
data = "".join(s1.split("\n")[1:])

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


# getting the count of a specific kmer,
# dividing by (length of sequence - length of kmer + 1)
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

# each box represents a k-mer
# get k-mers = 7
f7 = count_kmers(data, 7)
f7_prob = probabilities(f7, 7)
chaos_k7 = chaos_game_representation(f7_prob, 7)
chaos_k7

# run together for full plot
pylab.title('Chaos game representation for 7-mers')
k_mer_plot = pylab.imshow(chaos_k7, interpolation='nearest', cmap=cm.gray_r)
# pylab.show()
k_mer_plot

#_______________________________________________________________________________
# Calculating DCT
type(chaos_k7)
x = np.array(chaos_k7)
x
dct = dct(dct(x, type=2, norm='ortho'), type=3, norm='ortho')
dct

#_______________________________________________________________________________
from scipy.fftpack import dct, idct
import matplotlib.pyplot as plt
N = 100
t = np.linspace(0,20,N)
x = np.exp(-t/3)*np.cos(2*t)
y = dct(x, norm='ortho')
window = np.zeros(N)
window[:20] = 1
yr = idct(y*window, norm='ortho')
sum(abs(x-yr)**2) / sum(abs(x)**2)
#0.0010901402257
plt.plot(t, x, '-bx')
plt.plot(t, yr, 'ro')
window = np.zeros(N)
window[:15] = 1
yr = idct(y*window, norm='ortho')
sum(abs(x-yr)**2) / sum(abs(x)**2)
#0.0718818065008
plt.plot(t, yr, 'g+')
plt.legend(['x', '$x_{20}$', '$x_{15}$'])
plt.grid()
plt.show()

#_______________________________________________________________________________
# Finding Single Value Decomposition


# read image in grayscale
# img = cv2.imread('beach-2179624_960_720.jpg', 0)

# obtain svd
U, S, V = np.linalg.svd()

# inspect shapes of the matrices
print(U.shape, S.shape, V.shape)
