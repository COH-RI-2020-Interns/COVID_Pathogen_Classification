import numpy as np
import pandas as pd
import itertools
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.fft import fft, ifft
from random import sample

#Going to Test folders
folder_path = getcwd() + "/data3"

folders = sorted(listdir(folder_path))[0:17]

folder_dict = {}

# Going to Specific Virus Folders inside the Test folders
for folder in folders:
    folder_dict[folder] = listdir(f"{folder_path}/{folder}")


#Adding data to a JSON file
output_path = getcwd() + "/data3/JSON_Files"

with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}")

my_dict = json.load(f)
for test in my_dict:
    test = sorted(test)


#Dictionary of numerical representations
rep_dict = {"Int1":{"T":0,"t":0,"C":1,"c":1, "A":2,"a":2 ,"G":3, "g":3},
"Int2": {"T":1,"t":1,"C":2,"c":2, "A":3,"a":3 ,"G":4, "g":4},
"Real": {"T":-1.5,"t":-1.5,"C":0.5,"c":0.5, "A":1.5,"a":1.5 ,"G":-1.5, "g":-1.5},
"EIIP": {"T":0.1335,"t":0.1335,"C":0.1340,"c":0.1340, "A":0.1260,"a":0.1260 ,"G":0.0806, "g":0.0806},
"PP": {"T":1,"t":1,"C":1,"c":1, "A":-1,"a":-1 ,"G":-1, "g":-1},
"Paired Numeric": {"T":1,"t":1,"C":-1,"c":-1, "A":1,"a":1 ,"G":-1, "g":-1},
"Just A": {"T":0,"t":0,"C":0,"c":0, "A":1,"a":1 ,"G":0, "g":0},
"Just C": {"T":0,"t":0,"C":1,"c":1, "A":0,"a":0 ,"G":0, "g":0},
"Just G": {"T":0,"t":0,"C":0,"c":0, "A":0,"a":0 ,"G":1, "g":1},
"Just T": {"T":1,"t":1,"C":0,"c":0, "A":0,"a":0 ,"G":0, "g":0}}

input_dict = {'Int1':0,
'Int2':1,
'Real':2,
'EIIP':3,
'PP':4,
"Paired Numeric":5,
'Just A':6,
'Just C':7,
'Just G':8,
'Just T':9}
#____________________________________________________________________________________________________
# FUNCTIONS
# Finding the Average Magnitude of the Sequence using list comprehensions

def magnitude_avg(sequence):
    base_list = ["D", "K", "M", "N", "R", "S", "W", "Y", "H", "B","V"]
    for i in base_list:
        sequence = sequence.replace(i, "")
    dict_of_bases = [rep_dict[rep] for rep in rep_dict]
    numeric = [[dict_of_bases[i][base] for base in sequence] for i in range(0, len(rep_dict))]
    mag_avg_list = [np.average(abs(fft(np.array(list)))) for list in numeric]
    return mag_avg_list


# Calculating Entropy
# def entropy(sequence):
#     counts = Counter(sequence)
#     props = {key: counts[key] / sum(counts.values()) for key in counts}
#     products = {key: props[key]*np.log(props[key]) for key in props}
#     return -1 * sum(products.values())
#
# #Calculating Magtropy
# def magtropy(sequence):
#     list_magtropy = [avg/entropy(sequence) for avg in magnitude_avg(sequence)]
#     return list_magtropy

#sequence separation using list comprehension
def seq_separation_lst(sublevel, seq_num):
    file_path = getcwd()
    start_seq = [sample(list(SeqIO.parse((f"{file_path}/data3/{sublevel}/{file}"), "fasta")), len(list(SeqIO.parse((f"{file_path}/data3/{sublevel}/{file}"), "fasta")))) for file in my_dict[sublevel]]
    final_seq = [["".join([char for char in sequence.seq]) for sequence in list[0:seq_num]] for list in start_seq]
    seq_dict = {file_name[:-6]:list for (file_name,list) in zip(my_dict[sublevel], final_seq)}
    return seq_dict

# Saving Magtropy values to dictionary for specific sublevel using list comprehensions
def magtropy_dict(sublevel_dict):
    magtropy_values = [[(sublevel, magnitude_avg(sequence)[idx]) for sequence in sublevel_dict[sublevel]] for sublevel in sublevel_dict.keys()] #[0] retrieves the value within the list for each sequence
    taxonomic_level = pd.DataFrame((list(itertools.chain.from_iterable(magtropy_values))), columns = ["Sublevel Name", num_rep])#, "Int2", "Real", "EIIP", "PP", "Paired Numeric", "Just A", "Just C", "Just G", "Just T"])
    return taxonomic_level

#____________________________________________________________________________________________________

# DATA
taxonomic_level = input("Taxonomic level: ")
sublevel = seq_separation_lst(taxonomic_level, 100)
df_lst = []

# Preparing training data for supervised machine learning
num_rep = input("Numerical Rep [Int1, Int2, Real, EIIP, PP, Paired Numeric, Just A, Just C, Just G, Just T]: ")
idx = input_dict[num_rep]

df = magtropy_dict(sublevel)
df = df.drop(columns='Sublevel Name') # don't run for Int1
df


df_lst.append(df)
len(df_lst)

#only run after all numerical reps are added to df_lst
sublevel_df = pd.concat(df_lst, axis=1)

sublevel_df

sublevel_df['Sublevel Name'].value_counts()



saved_path = getcwd() + f"/Work_Remaining/num_rep_data_csvs"

sublevel_df.to_csv(saved_path + f"/num_reps_{taxonomic_level[2:]}_100.csv")
