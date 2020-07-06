import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from collections import Counter
import matplotlib.pyplot as plt

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


# Calculating Entropy
def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())


# Saving Entropy values to dictionary
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

entropy_dict['Test1']['Riboviria']

test_dict = {}
for test in entropy_dict:
    family_dict = {}
    for family in entropy_dict[test]:
        entropy_nums = []
        for value in entropy_dict[test][family]:
            entropy_nums.append(value[1])
            family_dict[family] = entropy_nums
    test_dict[test] = family_dict

len(test_dict['Test1'].keys())

for test in test_dict:
    for fam in test_dict[test].keys():
        entropy = test_dict[test][fam]
        plt.plot(entropy)
        plt.title("Entropy Values for " + fam)
        plt.legend(fam)

entropy1 = test_dict['Test1']['Riboviria']
len(entropy1)
entropy1 = entropy1[0:144]
len(entropy1)
entropy2 = test_dict['Test1']['Polyomaviridae']
len(entropy2)
plt.plot(entropy1)
plt.plot(entropy2)
plt.title("Entropy Values of Polyomaviridae & Riboviria")
plt.legend(["Riboviria", "Polyomaviridae"])
