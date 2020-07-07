import numpy as np
import pandas as pd
from Bio import SeqIO
import rapidjson
from os import getcwd, listdir

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


my_fasta['Test3a'].keys()
