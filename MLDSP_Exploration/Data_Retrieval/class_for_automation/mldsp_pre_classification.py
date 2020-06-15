import pandas as pd
from collections import Counter
from os import getcwd, listdir
from Bio import SeqIO
from scipy.fft import fft, ifft
from scipy import stats

class mlsdp_pre_classfication:
     def __init__(file_path1, file_path2):
         self.file_path = file_path
         self.file_path2 = file_path2

     # Extracting string from file
     def make_sequence1(number1):
         file_list = sorted(listdir(file_path))
         ribo_example = list(SeqIO.parse(f"{file_path}/{file_list[number1]}", "fasta"))
         count = len(ribo_example[0].seq)
         seq = "".join([char for char in ribo_example[0].seq])
         return "This is the first sequence: /n" + seq

     def make_sequence2(number2):
         file_list2 = sorted(listdir(file_path2))
         ribo_example = list(SeqIO.parse(f"{file_path2}/{file_list2[number2]}", "fasta"))
         count = len(ribo_example[0].seq)
         seq = "".join([char for char in ribo_example[0].seq])
         return "This is the second sequence: /n" + seq





# Getting 1st user input
test = input("Please input the 1st test you would like data from")
folder = input("Please input the 1st folder you would like data from")
folder
# Getting 2nd user input
test2 = input("Please input the 2nd test you would like data from")
folder2 = input("Please input the 2nd folder you would like data from")

# Getting file paths for both files
file_path = getcwd() + f"/data/{test}/{folder}"
file_path2 = getcwd() + f"/data/{test2}/{folder2}"
