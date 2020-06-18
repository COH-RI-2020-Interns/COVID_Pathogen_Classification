import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from itertools import permutations

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

list_1  = [1]
list_2 = [3]
list_1 + list_2


folder_dict.keys()
folder_dict.values()

# Adding data to a JSON file
output_path = getcwd() + "/data/JSON_Files"


with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = json.load(f)
my_dict.values()
my_dict["Test1"].keys()
# Getting all possible combinations of 2 for the fasta files
file_tuple_list = []

for folder in folders:
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        print(sub_folder)
        for file in listdir(f"{folder_path}/{folder}/{sub_folder}"):
            file_tuple_list.append((folder,sub_folder,file))
file_tuple_list




file_combos = list(permutations(file_tuple_list, 2))
file_combos



new_dict = {}
#each of the keys shows tuples, first element is the virus, second is file



for i in my_dict.keys():
    file_list = []
    for j in my_dict[i].keys():
        for file in my_dict[i][j]:
            file_list.append((j,file))
        new_dict[i] = file_list

new_dict

new_dict_2 = {}
for key in new_dict.keys():
    seq_perm = list(permutations(new_dict[key], 2))
    new_dict_2[key] =  seq_perm


new_dict_3 = {}
for key in new_dict_2.keys():
    for i,j in new_dict_2[key]:
        print(i)

len(my_dict)
file_tuple_list_practice = {}

for i in my_dict:
   file_tuple_list_practice  = {}
   file_tuple_list_practice = my_dict[i]
   print(file_tuple_list_practice)

file_perms

file_tuple_list


file_combos[6000000]


 #taking a while to go through everything
