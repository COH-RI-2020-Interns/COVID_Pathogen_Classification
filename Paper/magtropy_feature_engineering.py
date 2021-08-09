import numpy as np
import pandas as pd
import itertools
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.fft import fft, ifft
from random import sample




#Dictionary of numerical representations
rep_dict = {"Int1":{"T":0,"t":0,"C":1,"c":1, "A":2,"a":2 ,"G":3, "g":3},
"Int2": {"T":1,"t":1,"C":2,"c":2, "A":3,"a":3 ,"G":4, "g":4},
"Real": {"T":-1.5,"t":-1.5,"C":0.5,"c":0.5, "A":1.5,"a":1.5 ,"G":-1.5, "g":-1.5},
"EIIP": {"T":0.1335,"t":0.1335,"C":0.1340,"c":0.1340, "A":0.1260,"a":0.1260 ,"G":0.0806, "g":0.0806},
"PP": {"T":1,"t":1,"C":1,"c":1, "A":-1,"a":-1 ,"G":-1, "g":-1},
"Paired Numeric": {"T":1,"t":1,"C":-1,"c":-1, "A":1,"a":1 ,"G":-1, "g":-1},
"Just A": {"T":0,"t":0,"C":0,"c":0, "A":1,"a":1 ,"G":0, "g":0},
"Just C": {"T":0,"t":0,"C":1,"c":1, "A":0,"a":0 ,"G":0, "g":0},
"Just G": {"T":0,"t":0,"C":0,"c":0, "A":0,"a":0 ,"G":1, "g":1},
"Just T": {"T":1,"t":1,"C":0,"c":0, "A":0,"a":0 ,"G":0, "g":0}}


#____________________________________________________________________________________________________
# FUNCTIONS
# Finding the Average Magnitude of the Sequence using list comprehensions

def magnitude_avg(sequence):
    base_list = ["D", "K", "M", "N", "R", "S", "W", "Y", "H", "B","V"]
    for i in base_list:
        sequence = sequence.replace(i, "")
    dict_of_bases = [rep_dict[rep] for rep in rep_dict]
    numeric = [[dict_of_bases[i][base] for base in sequence] for i in range(0, len(rep_dict))]
    mag_avg_list = [np.average(abs(fft(np.array(list)))) for list in numeric]
    return mag_avg_list


# Calculating Entropy
def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())

#Calculating Magtropy
def mag_entropy(sequence):
    list_magtropy = [avg/entropy(sequence) for avg in magnitude_avg(sequence)]
    return list_magtropy

def entropy_mag(sequence):
    list_magtropy = [entropy(sequence)/avg for avg in magnitude_avg(sequence)]
    return list_magtropy
#____________________________________________________________________________________________________
