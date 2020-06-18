import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from itertools import permutations
from scipy.fft import fft, ifft
from scipy import stats



#Trying to get in JSON final_permutations
output_path = getcwd() + "/data/JSON_Files"
#with open(f"{output_path}/final_permutations.json", "r") as my_file:
    #json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[1]}", )

my_dict = json.load(f)
#JSON not loading
#_______________________________________________________________________________

#Using dictionary instead
folder_path = getcwd() + "/data"

folders = sorted(listdir(folder_path))[1:8]

folder_dict = {}

# Going to Specific Virus Folders inside the Test folders
for folder in folders:
    subfolder_dict = {}
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        subfolder_dict[sub_folder] = listdir(f"{folder_path}/{folder}/{sub_folder}")
    folder_dict[folder] = subfolder_dict

output_path = getcwd() + "/data/JSON_Files"


with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = json.load(f)


# Getting all possible combinations of 2 for the fasta files
file_tuple_list = []

for folder in folders:
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        for file in listdir(f"{folder_path}/{folder}/{sub_folder}"):
            file_tuple_list.append((folder,sub_folder,file))


file_combos = list(permutations(file_tuple_list, 2))

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
    seq_perm = list(permutations(new_dict[key], 2))
    new_dict_2[key] =  seq_perm


new_dict_2

new_dict_3 = {}
for key in new_dict_2.keys():
    file_list_2 = []
    for i,j in new_dict_2[key]:
        if(i[0] != j[0]):
            file_list_2.append((i,j))
    new_dict_3[key] = file_list_2

new_dict_3.keys()

new_dict_3["Test1"][0]


def make_sequence(path_of_file):
    start_seq = list(SeqIO.parse(f"{path_of_file}", "fasta"))
    count = len(start_seq[0].seq)
    final_seq = "".join([char for char in start_seq[0].seq])
    return final_seq






pearsons_dict = {}
for test in new_dict_3:
    for file1,file2 in new_dict_3[test]:
        file_path = getcwd() + f"/data/Test6/{file1[0]}/{file1[1]}"
        file_path2 = getcwd() + f"/data/Test6/{file2[0]}/{file2[1]}"
        seq1  = make_sequence(file_path)
        seq2 = make_sequence(file_path2)






pearsons_dict[test] = (file1, seq1, file2, seq2)













# Extracting string from file






polyomaviridae = make_sequence(file_list[1], file_path)
print(polyomaviridae)
riboviria = make_sequence(file_list2[1], file_path2)
print(riboviria)
len(riboviria)

# Length normalization
riboviria = riboviria[0:(len(polyomaviridae))]
len(riboviria)
riboviria
len(polyomaviridae)
polyomaviridae
# If we cut one sequence we could be taking away key genetic info

# Converting sequence based on JustA representation
dict_of_bases = {"T":1, "C":1, "A":-1, "G":-1}

def numerical(dna_strand):
    numeric = []
    for base in dna_strand:
        numeric.append(dict_of_bases[base])
    return np.array(numeric)

polyomaviridae_nums = np.array(numerical(polyomaviridae))
len(polyomaviridae_nums)
polyomaviridae_nums

riboviria_nums = np.array(numerical(riboviria))
len(riboviria_nums)
riboviria_nums


#Calculating discrete numerical representation
𝐹𝑖(𝑘)=∑𝑗=0𝑝−1𝑓(𝑆𝑖(𝑗))⋅𝑒(−2𝜋𝑖/𝑝)𝑘𝑗
#Results  = [6, -1 -3i , 0, 1+3i]
#Logic =
#to get 6 do (1) + (2) + (3)
#to get -1 -3i do 1 + 3*e^(-pi(i)/2) + 2*e^(-pi(i)) + 0
#to get 0 do 1 + 3*e^(-pi(i)) + 2*e^(-2pi(i)) + 0
#to get -1+3i 1 + 3*e^(-3/2(pi(i))) + 2*e^(-3(pi)(i)) + 0

dft_polyomaviridae = fft(polyomaviridae_nums)
len(dft_polyomaviridae)
dft_polyomaviridae
dft_riboviria = fft(riboviria_nums)
len(dft_riboviria)
dft_riboviria

#Finding magnitude
mag_list = []
def mag(list_of_dtf):
    mag_list = abs(list_of_dtf)
    return mag_list


mag_polyomaviridae = mag(dft_polyomaviridae)
mag_polyomaviridae
len(mag_polyomaviridae)

mag_riboviria = mag(dft_riboviria)
mag_riboviria
len(mag_riboviria)

#DTF magnitude spectra graph PP representation
poly_50 = mag_polyomaviridae[0:50]
ribo_50 = mag_riboviria[0:50]
plt.plot(poly_50)
plt.plot(ribo_50)
plt.title("First 50 DFT Magnitude Spectra PP representation")
plt.legend(["Poly_50", "Ribo_50"])


#Discrete digital signal graph PP representation
polyomaviridae_nums_50 = polyomaviridae_nums[0:50]
riboviria_nums_50 = riboviria_nums[0:50]
plt.plot(polyomaviridae_nums_50)
plt.plot(riboviria_nums_50)
plt.title("First 50 discrete digital signal PP representation")
plt.legend(["Poly_dds_50", "Ribo_dds_50"])
#Purines are -1 (AG) and Pyrimidines are 1 (TC)


#Finding Pearson's Correlation Coefficient
stats.pearsonr(mag_polyomaviridae, mag_riboviria)
# PCC is extremely sensitive to extreme values
# The p-value measures the percent risk that the data seems to correlate when there is actually no correlation
# testing if caused by chance
# For p-values equal to or below 0.05, the data is stastically significant
# For p-values above 0.05, the data is not stastically significant and should not be considered
