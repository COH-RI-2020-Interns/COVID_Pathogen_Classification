import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from itertools import combinations

#fake_files = [f"file_{i+1}.fasta" for i in range(100)]

file_path = getcwd() + "/data/Test1/Riboviria"

#for file_name in fake_files:
    #system(f"touch {file_path}/{file_name}")

folder_path = getcwd() + "/data"

folders = sorted(listdir(folder_path))[2:9]
folders

folder_dict = {}

for folder in folders:
    subfolder_dict = {}
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        subfolder_dict[sub_folder] = listdir(f"{folder_path}/{folder}/{sub_folder}")
    folder_dict[folder] = subfolder_dict

folder_dict.keys()

output_path = getcwd() + "/data/JSON_Files"

folder_json = json.dumps(folder_dict)

#system(f"touch {output_path}/fasta_files.json")

with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)


coords = [("ribo", 4,5), ("sarbecovirus", 2,8), ("ribo", 1,9), ("sarbecovirus", 20,8), ("ribo", 7,19)]

coord_pairs = combinations(coords, 2)

list(coord_pairs)[0]
