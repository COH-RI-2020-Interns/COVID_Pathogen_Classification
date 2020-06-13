import numpy as np
from collections import Counter
from os import getcwd, listdir


from Bio import SeqIO

# Loading in the Data
file_path = getcwd() + "/Shannon entropy data"
file_path
file = listdir(file_path)[1]
file
ribo_example = list(SeqIO.parse(f"{file_path}/{file}", "fasta"))

len(ribo_example[0].seq)

seq = "".join([char for char in ribo_example[0].seq])

def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())

entropy(seq)
