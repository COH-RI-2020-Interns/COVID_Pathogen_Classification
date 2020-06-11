import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from os import getcwd, listdir
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import Counter

file_path = getcwd() + "/COVID_Pathogen_Classification/data/Test1/Riboviria"

# Reading in the Data as a List Comprehension
ribovirus_example =  [line.replace("\n", "") for line in open(f"{file_path}/{listdir(file_path)[0]}", "r").readlines()]

# Transforming Each Sequence into a BioPython Seq Object
ribovirus_object = [Seq(line, generic_dna) for line in ribovirus_example]

# Getting Frequency of Each Base in DNA
dna_base_frequencies = [Counter(seq) for seq in ribovirus_example]
dna_base_frequencies
totals = [sum(list(freq.values())) for freq in dna_base_frequencies]

dna_base_frequencies[0]

freq_prop = [{key : freq[key] / sum(list(freq.values())) for key in freq.keys()} for freq in dna_base_frequencies]

entropy = [-1 * sum({key : freq[key]*np.log(freq[key]) for key in freq.keys()}.values()) for freq in freq_prop]
entropy
