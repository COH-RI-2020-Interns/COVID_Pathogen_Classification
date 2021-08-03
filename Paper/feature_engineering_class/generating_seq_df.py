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



#
# def seq_separation_dict(taxonomic_level):
#     file_path = getcwd()
#     start_seq = [sample(list(SeqIO.parse((f"{file_path}/data_final/{taxonomic_level}/{file}"), "fasta")), len(list(SeqIO.parse((f"{file_path}/data_final/{taxonomic_level}/{file}"), "fasta")))) for file in my_dict[taxonomic_level]]
#     final_seq = [["".join([char for char in sequence.seq]) for sequence in list] for list in start_seq]
#     seq_dict = {file_name[:-6]:list for (file_name,list) in zip(my_dict[taxonomic_level], final_seq)}
#     return seq_dict


tax_lvl_list = ['0_Covid', '1_Realm', '2_Kingdom', '3_Phylum', '4_Class', '5_Order', '6_Suborder', '8_Genus', '9_Subgenus']
sequence_dict = {}
# for tax_lvl in tax_lvl_list:
#     sequence_dict.update(seq_separation_dict(tax_lvl))


def seq_separation_dict_ctrl(sublevel, seq_num):
    file_path = getcwd()
    start_seq = [sample(list(SeqIO.parse((f"{file_path}/data_final/{sublevel}/{file}"), "fasta")), len(list(SeqIO.parse((f"{file_path}/data_final/{sublevel}/{file}"), "fasta")))) for file in my_dict[sublevel]]
    final_seq = [["".join([char for char in sequence.seq]) for sequence in list[0:seq_num]] for list in start_seq]
    seq_dict = {file_name[:-6]:list for (file_name,list) in zip(my_dict[sublevel], final_seq)}
    return pd.DataFrame.from_dict(seq_dict)

seq_separation_dict_ctrl("6_Suborder", 1)
























# seq_df = pd.DataFrame
# for tax_lvl in tax_lvl_list:
#     seq_df[tax_lvl[2:]] = seq_separation_dict_ctrl(tax_lvl, 1)


# for tax_lvl in tax_lvl_list:
#     sequence_dict.update(seq_separation_dict_ctrl(tax_lvl, 1))









# taxonomic_level = input("Taxonomic level: ")
# sublevel = seq_separation_lst(taxonomic_level, 100)
#
# sublevel_df = magtropy_dict(sublevel)
#
# sublevel_df['Sublevel_Name'].value_counts()
#
# sublevel_df
