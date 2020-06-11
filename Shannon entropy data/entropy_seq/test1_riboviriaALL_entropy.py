import numpy as np
from collections import Counter
from os import getcwd, listdir


from Bio import SeqIO

# Loading in the Data
file_path = getcwd() + "/COVID_Pathogen_Classification/Shannon entropy data/data/Test1/Riboviria"
file_path

for i in range(len(file_path)-1):
    file = listdir(file_path)[i]

    for f in range(len(file) - 1):
        ribo_example = list(SeqIO.parse(f"{file_path}/{file}", "fasta"))
        #ribo_example
        count = len(ribo_example[0].seq)
        #print(count)
        seq = "".join([char for char in ribo_example[0].seq])
        #seq
    def entropy(sequence):
        counts = Counter(sequence)
        props = {key: counts[key] / sum(counts.values()) for key in counts
        products = {key: props[key]*np.log(props[key]) for key in props}
        return -1 * sum(products.values())

        entropy(seq)
