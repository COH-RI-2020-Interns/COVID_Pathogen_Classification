import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from collections import Counter
from Bio import SeqIO


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

entropy_values = []

my_dict.keys()
my_dict['Test1'].keys()
my_dict["Test1"]["Polyomaviridae"]




my_dict
my_dict[]
file_path_1 = getcwd() + "/data"
for folder in my_dict.keys():
    for sub_folder in my_dict[folder].keys():
        #print(sub_folder)
        for file in my_dict[folder][sub_folder]:
            file_path_1 = listdir(f"/data/{folder}/{sub_folder}/{file}")
            ribo_example = list(SeqIO.parse(f"{my_dict[folder][sub_folder][file]}", "fasta"))
            print(ribo_example)
            count = len(ribo_example[0].seq)
            seq = "".join([char for char in ribo_example[0].seq])
            entropy_values.append(entropy(seq))

entropy_values #more entropy = more info (ML) = more uncertainty

def Avg(lst):
    return sum(lst) / len(lst)

Avg(entropy_values)
