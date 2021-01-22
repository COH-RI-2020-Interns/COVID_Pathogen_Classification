import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from itertools import combinations

fake_files = [f"file_{i+1}.fasta" for i in range(100)]

file_path = getcwd() + "/data/Fake_Data/Fake_Virus_1"

for file_name in fake_files:
    system(f"touch {file_path}/{file_name}")

fake_folder_path = getcwd() + "/data"

fake_folders = sorted(listdir(fake_folder_path))[1:3]

fake_folder_dict = {}

for folder in fake_folders:
    fake_subfolder_dict = {}
    for sub_folder in listdir(f"{fake_folder_path}/{folder}"):
        fake_subfolder_dict[sub_folder] = listdir(f"{fake_folder_path}/{folder}/{sub_folder}")
    fake_folder_dict[folder] = fake_subfolder_dict

fake_folder_dict.keys()

output_path = getcwd() + "/data/JSON_Files"

folder_json = json.dumps(fake_folder_dict)

system(f"touch {output_path}/fasta_files.json")

with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(fake_folder_dict, my_file)


coords = [("ribo", 4,5), ("sarbecovirus", 2,8), ("ribo", 1,9), ("sarbecovirus", 20,8), ("ribo", 7,19)]

coord_pairs = combinations(coords, 2)

list(coord_pairs)[0]
