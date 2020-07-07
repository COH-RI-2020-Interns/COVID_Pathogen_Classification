import numpy as np
import pandas as pd
from Bio import SeqIO
import rapidjson
from os import getcwd, listdir
from itertools import permutations, combinations
from scipy.fft import fft, ifft
from scipy import stats
import pywt

data_path = getcwd() + "/data/JSON_Files"

with open(f"{data_path}/{listdir(data_path)[0]}", "r") as my_file:
    my_fasta = rapidjson.load(my_file)

my_fasta.keys()

cluster_info = []
cluster_name = []
for key in my_fasta['Test3a'].keys():
    for file_name in my_fasta['Test3a'][key]:
         cluster_name.append(key)
         cluster_info.append(len(my_fasta['Test3a'][key]))




#Using dictionary instead
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

output_path = getcwd() + "/data/JSON_Files"


with open(f"{output_path}/fasta_files.json", "w") as my_file:
    rapidjson.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = rapidjson.load(f)


# Getting all possible combinations of 2 for the fasta files
file_tuple_list = []

for folder in folders:
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        for file in listdir(f"{folder_path}/{folder}/{sub_folder}"):
            file_tuple_list.append((sub_folder,file))


file_combos = list(combinations(file_tuple_list, 2))


new_dict = {}
#each of the keys shows tuples, first element is the virus, second is file

for i in my_dict.keys():
    file_list = []
    for j in my_dict[i].keys():
        for file in my_dict[i][j]:
            file_list.append((j,file))
        new_dict[i] = file_list



new_dict_2 = {}
for key in new_dict.keys():
    seq_perm = list(combinations(new_dict[key], 2))
    new_dict_2[key] =  seq_perm


len(new_dict_2["Test3a"])
82 * 82
