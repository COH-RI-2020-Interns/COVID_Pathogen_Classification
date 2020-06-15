import numpy as np
from collections import Counter
from os import getcwd, listdir

#riboviria all
from Bio import SeqIO

# Loading in the Data
file_path = getcwd() + "/data/Test1/Riboviria"
listdir(file_path)

def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())

entropy_values = []

#for i in range(len(listdir(file_path))-1):
    #file = listdir(file_path)[i]
    #file


for file in listdir(file_path):
    ribo_example = list(SeqIO.parse(f"{file_path}/{file}", "fasta"))

    count = len(ribo_example[0].seq)

    seq = "".join([char for char in ribo_example[0].seq])

        #print(entropy(seq))
    entropy_values.append(entropy(seq))

entropy_values #more entropy = more info (ML) = more uncertainty

def Avg(lst):
    return sum(lst) / len(lst)

Avg(entropy_values)
