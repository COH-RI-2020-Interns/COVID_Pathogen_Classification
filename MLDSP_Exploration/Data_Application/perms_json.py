import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from itertools import permutations
from scipy.fft import fft, ifft
from scipy import stats
import pywt



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

folders = sorted(listdir(folder_path))[2:9]
folders
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

new_dict_3 = {}
for key in new_dict_2.keys():
    file_list_2 = []
    for i,j in new_dict_2[key]:
        if(i[0] != j[0]):
            file_list_2.append((i,j))
    new_dict_3[key] = file_list_2

#___________________________________________________________________________
#Using the dictionary to apply to our code

def make_sequence(path_of_file):
    start_seq = list(SeqIO.parse(f"{path_of_file}", "fasta"))
    count = len(start_seq[0].seq)
    final_seq = "".join([char for char in start_seq[0].seq])
    return final_seq

dict_of_bases = {"T":1, "C":1, "A":-1, "G":-1}

def numerical_pp(dna_strand):
    numeric = []
    for base in dna_strand:
        numeric.append(dict_of_bases[base])
    return np.array(numeric)


def normalization(numerical1, numerical2):
    if len(numerical1)>len(numerical2):
        dna_seq = np.array(numerical2)
        numer = len(numerical1)- len(numerical2)
        if(numer%2 != 0):
            numer = numer - 0.5
            pad_width = numer/2
            numerical2 = pywt.pad(dna_seq, pad_width, "antisymmetric")
            numerical1 = numerical1[0:len(numerical1)-1]
        else:
            pad_width = numer/2
            numerical2 = pywt.pad(dna_seq,pad_width, "antisymmetric")
    else:
        dna_seq = np.array(numerical1)
        numer = len(numerical2)- len(numerical1)
        if(numer%2 != 0):
            numer = numer - 0.5
            pad_width = numer/2
            numerical1 = pywt.pad(dna_seq, pad_width, "antisymmetric")
            numerical2 = numerical2[0:len(numerical2)-1]
        else:
            pad_width = numer/2
            numerical1 = pywt.pad(dna_seq,pad_width, "antisymmetric")

    return numerical1,numerical2




pearsons_dict = {}
for test in new_dict_3:
    list_sequences= []
    for file1,file2 in new_dict_3[test]:
        file_path = getcwd() + f"/data/{test}/{file1[0]}/{file1[1]}"
        file_path2 = getcwd() + f"/data/{test}/{file2[0]}/{file2[1]}"
        seq1  = make_sequence(file_path)
        seq2 = make_sequence(file_path2)
        pp1 = numerical_pp(seq1)
        pp2 = numerical_pp(seq2)
        pp1_norm = normalization(pp1,pp2)[0]
        pp2_norm = normalization(pp1,pp2)[1]
        #print(len(pp1),len(pp1_norm), file1[1],len(pp2),len(pp2_norm),file2[1])
        fft_1 = fft(pp1_norm)
        fft_2 = fft(pp2_norm)
        mag_1 = abs(fft_1)
        mag_2 = abs(fft_2)
        pcc = stats.pearsonr(mag_1, mag_2)
        list_sequences.append((file1[1],file2[1],pcc))
    pearsons_dict[test] = list_sequences


#pearsons_dict












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
