from os import getcwd, listdir
import pandas as pd
import numpy as np
import json
from random import sample
from Bio import SeqIO

#Going to Test folders
folder_path = getcwd() + "/data_final"

folders = sorted(listdir(folder_path))[0:9]
folders
folder_dict = {}

# Going to Specific Virus Folders inside the Test folders
for folder in folders:
    folder_dict[folder] = listdir(f"{folder_path}/{folder}")


#Adding data to a JSON file
output_path = getcwd() + "/data_final/JSON_Files"

with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}")

my_dict = json.load(f)
for test in my_dict:
    test = sorted(test)



tax_lvl_list = ['0_Covid', '1_Realm', '2_Kingdom', '3_Phylum', '4_Class', '5_Order', '6_Suborder', '8_Genus', '9_Subgenus']
sequence_dict = {}

taxonomic_level = input("taxonomic level: ")

def seq_separation_dict_ctrl(taxonomic_level):
    file_path = getcwd()
    start_seq = [sample(list(SeqIO.parse((f"{file_path}/data_final/{taxonomic_level}/{file}"), "fasta")), len(list(SeqIO.parse((f"{file_path}/data_final/{taxonomic_level}/{file}"), "fasta")))) for file in my_dict[taxonomic_level]]
    final_seq = [["".join([char for char in sequence.seq]) for sequence in list] for list in start_seq]
    seq_dict = {file_name[:-6]:list for (file_name,list) in zip(my_dict[taxonomic_level], final_seq)}
    return pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in seq_dict.items() ]))


sublevel_df = seq_separation_dict_ctrl(taxonomic_level)
sublevel_df

saved_path = getcwd() + f"/Paper/feature_engineering_class/training_seq_csvs"

sublevel_df.to_csv(saved_path + f"/{taxonomic_level[2:]}_seqs.csv")



realm_df = pd.read_csv(f"{saved_path}/Realm_seqs.csv")
kingdom_df = pd.read_csv(f"{saved_path}/Kingdom_seqs.csv")
phylum_df = pd.read_csv(f"{saved_path}/Phylum_seqs.csv")
class_df = pd.read_csv(f"{saved_path}/Class_seqs.csv")
order_df = pd.read_csv(f"{saved_path}/Order_seqs.csv")
suborder_df = pd.read_csv(f"{saved_path}/Suborder_seqs.csv")
genus_df = pd.read_csv(f"{saved_path}/Genus_seqs.csv")
subgenus_df = pd.read_csv(f"{saved_path}/Subgenus_seqs.csv")



training_seqs_df = realm_df.merge(kingdom_df.merge(phylum_df.merge(class_df.merge(order_df.merge(suborder_df.merge(genus_df.merge(subgenus_df, how='left'), how='left'), how='left'), how='left'), how='left'), how='left'), how='left')


training_seqs_df.to_csv(saved_path + f"/master_seqs_df.csv")


saved_path2 = getcwd() + f"/Paper/feature_engineering_class/testing_seq_csv"



csvs = listdir(saved_path)
csvs

max = 0
val = csvs[0]
for i in csvs:
    new = len(pd.read_csv(f"{saved_path}/{i}"))
    if new > max:
        max = new
        val = i

print(val, max)


def norm_df_length():
    length = 69494
    for file in csvs:
        df = pd.read_csv(f"{saved_path}/{file}")
        df.drop(column = 'Unamed: 0')
        new_len = length - len(file)

list_lens =[]
for i in csvs:
    list_lens.append((len(pd.read_csv(f"{saved_path}/{i}")), i))
    list_lens.sort()

print(list_lens)
