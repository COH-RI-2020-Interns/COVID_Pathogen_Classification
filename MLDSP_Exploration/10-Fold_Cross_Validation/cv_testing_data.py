import numpy as np
import pandas as pd
import json
from sklearn.model_selection import train_test_split
from sklearn import datasets
from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn import metrics
from os import getcwd, listdir, system
from itertools import permutations
from scipy.fft import fft, ifft
from scipy import stats
from Bio import SeqIO
import pywt



#Using dictionary instead
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

output_path = getcwd() + "/data/JSON_Files"


with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = json.load(f)
my_dict
______________________________________________________________________
#Using the dictionary to apply to our code

def make_sequence(path_of_file):
    start_seq = list(SeqIO.parse(f"{path_of_file}", "fasta"))
    count = len(start_seq[0].seq)
    final_seq = "".join([char for char in start_seq[0].seq])
    return final_seq

dict_of_bases = {"T":1,"t":1,"C":1,"c":1, "A":-1,"a":-1 ,"G":-1, "g":-1}

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



#___________________________________________________________________
# Printing PCC data by Test


new_dict_4 = {}
new_dict_5 = {}
for folder in my_dict['Test3b']:
    list_sequences = []
    for file in my_dict["Test3b"][folder]:
        file_path = getcwd() + f"/data/Test3b/{folder}/{file}"
        seq  = make_sequence(file_path)
        pp = numerical_pp(seq)
        fft = np.fft.fft(pp)
        mag = abs(fft)
        list_sequences.append(mag)
    new_dict_4[folder] = list_sequences


new_dict_4["Deltacoronavirus"]
new_dict_4["Alphacoronavirus"]
new_dict_4["Betacoronavirus"]



#___________________________________________________________________

#train_test_split
#Using delta and alphacoronavirus data from Test3b
delta = new_dict_4["Deltacoronavirus"]
alpha = new_dict_4["Alphacoronavirus"]
test3_data = [delta]

#making an array of all the names for the first column (Delta/Alpha)
name1 = []
for i in range(len(delta)):
    name1.append("Deltacoronavirus")
for i in range(len(alpha)):
    name1.append("Alphacoronavirus")

#dataframe has all the values of the magnitudes split up into base pairs
#inserting family name
df = pd.DataFrame(data= delta + alpha)
df.insert(0, "Family", name1)

#Setting X to be the family name, y to be the magnitudes 
X = df["Family"]
y = df.drop(columns = ["Family"])


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)

clf = svm.SVC(kernel='linear', C=1).fit(X_train, y_train)
clf.score(X_test, y_test)




#essentially we are trying to find how correlated the 2 sequences are
#
