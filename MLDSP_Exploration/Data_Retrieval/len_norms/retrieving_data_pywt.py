import numpy as np
import pandas as pd
from collections import Counter
from os import getcwd, listdir
from Bio import SeqIO
from scipy.fft import fft, ifft
from scipy import stats
import pywt

# Getting 1st user input
test = input("Please input the 1st test you would like data from")
folder = input("Please input the 1st folder you would like data from")

# Getting 2nd user input
test2 = input("Please input the 2nd test you would like data from")
folder2 = input("Please input the 2nd folder you would like data from")

# Getting file paths for both files
file_path = getcwd() + f"/data/{test}/{folder}"
file_path2 = getcwd() + f"/data/{test2}/{folder2}"

len(listdir(file_path))
len(listdir(file_path2))

# Sorting folders and selecting file
file_list = sorted(listdir(file_path))
file_list[0]
file_list2 = sorted(listdir(file_path2))
file_list2[0]
type(file_list)

# Extracting string from file
def make_sequence(file, path_of_file):
    ribo_example = list(SeqIO.parse(f"{path_of_file}/{file}", "fasta"))
    count = len(ribo_example[0].seq)
    seq = "".join([char for char in ribo_example[0].seq])
    return seq

polyomaviridae = make_sequence(file_list[0], file_path)
print(polyomaviridae)
riboviria = make_sequence(file_list2[0], file_path2)
print(riboviria)
len(riboviria)

# Length normalization
len(riboviria)
len(polyomaviridae)
type(polyomaviridae)
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
#run strand 1 and then rerun for strand 2



if len(riboviria_nums)>len(polyomaviridae_nums):
    dna_seq = np.array(polyomaviridae_nums)
    numer = len(riboviria)- len(polyomaviridae)
    if(numer%2 != 0):
        numer = numer - 0.5
        pad_width = numer/2
        polyomaviridae_nums = pywt.pad(dna_seq, pad_width, "antisymmetric")
        riboviria_nums = riboviria_nums[0:len(riboviria_nums)-1]
    else:
        pad_width = numer/2
        polyomaviridae_nums = pywt.pad(dna_seq,pad_width, "antisymmetric")
else:
    dna_seq = np.array(riboviria_nums)
    numer = len(polyomaviridae_nums)- len(riboviria_nums)
    if(numer%2 != 0):
        numer = numer - 0.5
        pad_width = numer/2
        riboviria_nums = pywt.pad(dna_seq, pad_width, "antisymmetric")
        polyomaviridae_nums = numerical2[0:len(numerical2)-1]
    else:
        pad_width = numer/2
        riboviria_nums = pywt.pad(dna_seq,pad_width, "antisymmetric")



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


#Finding Pearson's Correlation Coefficient
stats.pearsonr(mag_polyomaviridae, mag_riboviria)
# PCC is extremely sensitive to extreme values
# The p-value measures the percent risk that the data seems to correlate when there is actually no correlation
# testing if caused by chance
# For p-values equal to or below 0.05, the data is stastically significant
# For p-values above 0.05, the data is not stastically significant and should not be considered
