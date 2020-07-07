import numpy as np
import pandas as pd
from Bio import SeqIO
import rapidjson
from os import getcwd, listdir

data_path = getcwd() + "/data/JSON_Files"

with open(f"{data_path}/{listdir(data_path)[1]}", "r") as perm:
    perm_vals = rapidjson.load(perm)

perm_vals.keys()

perm_vals['Test3a'][0]

[np.array((file1, file2)), np.array((file1, file3))]

pearson_corr = []
for perm in perms:
    pearson_corr.append(pearsons(perm))

pearson_corr = np.array(pearson_corr).reshape(82,82)
