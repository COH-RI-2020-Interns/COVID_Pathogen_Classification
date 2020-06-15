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



     # Sorting folders and selecting file
     file_list = sorted(listdir(file_path))
     file_list2 = sorted(listdir(file_path2))


     # Extracting string from file
     def make_sequence1(file_path, number1):
         ribo_example = list(SeqIO.parse(f"{file_path}/{file_list[number1]}", "fasta"))
         count = len(ribo_example[0].seq)
         seq = "".join([char for char in ribo_example[0].seq])
         return "This is the first sequence: /n" + seq

     def make_sequence2(file_path2, number2):
         ribo_example = list(SeqIO.parse(f"{file_path2}/{file_list2[number2]}", "fasta"))
         count = len(ribo_example[0].seq)
         seq = "".join([char for char in ribo_example[0].seq])
         return "This is the second sequence: /n" + seq



first_test = mlsdp_pre_classfication("Test1", "Polyomaviridae", "0", "Test1", "Riboviria", "0")
first_test.show_full_name()
