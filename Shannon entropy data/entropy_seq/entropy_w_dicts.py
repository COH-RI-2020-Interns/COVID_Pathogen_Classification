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

entropy_values = []

my_dict.keys()
my_dict['Test1'].keys()
my_dict["Test1"]["Polyomaviridae"]


file_path_1 = getcwd()
file_path_1
count_1 = 0
entropy_values = []
for folder in my_dict.keys():
    #print(folder)
    for sub_folder in my_dict[folder].keys():
        #print(sub_folder)
        for file in my_dict[folder][sub_folder]:
            count_1 = count_1 +1
            #file_path_1 = listdir(f"/data/{folder}/{sub_folder}/{file}")
            ribo_example = list(SeqIO.parse((f"{file_path_1}/data/{folder}/{sub_folder}/{file}"), "fasta"))
            #print(ribo_example)
            count = len(ribo_example[0].seq)
            seq = "".join([char for char in ribo_example[0].seq])
            entropy_values.append((folder, sub_folder, file, entropy(seq)))

entropy_values #more entropy = more info (ML) = more uncertainty
len(entropy_values)
count_1





# Getting all possible combinations of 2 for the fasta files
file_tuple_list = []

for folder in folders:
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        for file in listdir(f"{folder_path}/{folder}/{sub_folder}"):
            file_tuple_list.append((folder,sub_folder,file))

new_dict = {}
#each of the keys shows tuples, first element is the virus, second is file

for i in my_dict.keys():
    file_list = []
    for j in my_dict[i].keys():
        for file in my_dict[i][j]:
            file_list.append((j,file))
        new_dict[i] = file_list

new_dict['Test1']
len(new_dict['Test1'])

new_dict_2 = {}
for key in new_dict.keys():
    seq_perm = list(permutations(new_dict[key], 2))
    new_dict_2[key] =  seq_perm


new_dict_2

new_dict_3 = {}
for key in new_dict_2.keys():
    file_list_2 = []
    for i,j in new_dict_2[key]:
        if(i[0] != j[0]):
            file_list_2.append((i,j))
    new_dict_3[key] = file_list_2

new_dict_3
