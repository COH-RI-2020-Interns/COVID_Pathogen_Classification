import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from collections import Counter
from Bio import SeqIO


# Going to Test folders
folder_path = getcwd() + "/data"

folders = sorted(listdir(folder_path))[2:9]
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
my_dict['Test1'].keys()

# Getting all possible combinations of 2 for the fasta files
#file_tuple_list = []

#file_tuple_list
# Loading in the Data

def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())

my_dict.keys()
my_dict['Test1'].keys()
my_dict["Test1"]["Polyomaviridae"]


entropy_dict = {}
for test in my_dict.keys():
    entropy_values = []
    temp_entropy_dict = {}
    for family in  my_dict[test].keys():
        for file in my_dict[test][family]:
            start_seq = list(SeqIO.parse((f"{file_path_1}/data/{test}/{family}/{file}"), "fasta"))
            count = len(start_seq[0].seq)
            final_seq = "".join([char for char in start_seq[0].seq])
            entropy_values.append((file, entropy(seq)))
            temp_entropy_dict[family] = entropy_values
    entropy_dict[test] = temp_entropy_dict

temp_entropy_dict
