import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
import collections


# Going to Test folders
folder_path = getcwd() + "/data"

folders = sorted(listdir(folder_path))[1:8]
folders

folder_dict = {}

# Going to Specific Virus Folders inside the Test folders
for folder in folders:
    subfolder_dict = {}
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        subfolder_dict[sub_folder] = listdir(f"{folder_path}/{folder}/{sub_folder}")
    folder_dict[folder] = subfolder_dict

folder_dict.keys()

# Adding data to a JSON file
output_path = getcwd() + "/data/JSON_Files"


with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = json.load(f)


# Extracting Sequence from Files

def make_sequence(path_of_file):
    start_seq = list(SeqIO.parse(f"{path_of_file}", "fasta"))
    count = len(start_seq[0].seq)
    final_seq = "".join([char for char in start_seq[0].seq])
    return final_seq

# Counting total k-mers possible
def getKmers(sequence, k):
    return [sequence[x:x+k].lower() for x in range(len(sequence) - k + 1)]

seq = "GAGAGACAAAGTTCAAAGGGCTATACAACCCCTGAATAGTAACAAAATACAGAAAAACCATAAAATTATAAAAATAACTAATCTGATCATCTAAATTTGACTAATTGGAAATAGCCGAACTCTACGGAGATGTAGGCGTCCGAACTCCACGGAGACGTAGGACAAAATTCTGCCGAACCCCAGACCATCGGGGACGTAGGCGTCTAATTTGTTTTTTTAATATTTTAC"
len(seq)
def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())


# Find entropy of k-mers
def entropy_k(kmer):
    kmer_lst = []
    for sequence in kmer:
        counts = Counter(sequence)
        props = {key: counts[key] / sum(counts.values()) for key in counts}
        products = {key: props[key]*np.log(props[key]) for key in props}
        entropy_kmer = -1 * sum(products.values())
        kmer_lst.append(entropy_kmer)
    return np.average(kmer_lst)

entropy(k)
entropy_k(getKmers(seq, 5))

file_path_1 = getcwd()
entropy_dict = {}
for test in my_dict.keys():
    temp_entropy_dict = {}
    for family in  my_dict[test].keys():
        entropy_values = []
        for file in my_dict[test][family]:
            start_seq = list(SeqIO.parse((f"{file_path_1}/data/{test}/{family}/{file}"), "fasta"))
            count = len(start_seq[0].seq)
            final_seq = "".join([char for char in start_seq[0].seq])
            entropy_values.append((file, entropy(final_seq)))
            temp_entropy_dict[family] = entropy_values
    entropy_dict[test] = temp_entropy_dict



entropy_dict["Test1"]["Polyomaviridae"]
entropy_dict['Test1']['Riboviria']
