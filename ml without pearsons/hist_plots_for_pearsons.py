import numpy as np
import pandas as pd
import json
from Bio import SeqIO
from os import getcwd, listdir, system
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.fft import fft, ifft
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, classification_report

#Going to Test folders
folder_path = getcwd() + "/data"

folders = sorted(listdir(folder_path))[2:10]
folders

folder_dict = {}

# Going to Specific Virus Folders inside the Test folders
for folder in folders:
    subfolder_dict = {}
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        subfolder_dict[sub_folder] = listdir(f"{folder_path}/{folder}/{sub_folder}")
    folder_dict[folder] = subfolder_dict

folder_dict.keys()

#Adding data to a JSON file
output_path = getcwd() + "/data/JSON_Files"


with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = json.load(f)
for test in my_dict:
    test = sorted(test)

# Calculating Entropy
def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())

def magnitude_avg(sequence, representation = "PP"):
    if representation == "Int1":
        dict_of_bases = {"T":0,"t":0,"C":1,"c":1, "A":2,"a":2 ,"G":3, "g":3}
    elif representation == "Int2":
        dict_of_bases = {"T":1,"t":1,"C":2,"c":2, "A":3,"a":3 ,"G":4, "g":4}
    elif representation == "Real":
        dict_of_bases = {"T":-1.5,"t":-1.5,"C":0.5,"c":0.5, "A":1.5,"a":1.5 ,"G":-1.5, "g":-1.5}
    elif representation == "Atomic":
        dict_of_bases = {"T":6,"t":6,"C":58,"c":58, "A":70,"a":70 ,"G":78, "g":78}
    elif representation == "EIIP":
        dict_of_bases = {"T":0.1335,"t":0.1335,"C":0.1340,"c":0.1340, "A":0.1260,"a":0.1260 ,"G":0.0806, "g":0.0806}
    elif representation == "PP":
        dict_of_bases = {"T":1,"t":1,"C":1,"c":1, "A":-1,"a":-1 ,"G":-1, "g":-1}
    elif representation == "Paired Numeric":
        dict_of_bases = {"T":1,"t":1,"C":-1,"c":-1, "A":1,"a":1 ,"G":-1, "g":-1}
    elif representation == "Just A":
        dict_of_bases = {"T":0,"t":0,"C":0,"c":0, "A":1,"a":1 ,"G":0, "g":0}
    elif representation == "Just C":
        dict_of_bases = {"T":0,"t":0,"C":1,"c":1, "A":0,"a":0 ,"G":0, "g":0}
    elif representation == "Just G":
        dict_of_bases = {"T":0,"t":0,"C":0,"c":0, "A":0,"a":0 ,"G":1, "g":1}
    elif representation == "Just T":
        dict_of_bases = {"T":1,"t":1,"C":0,"c":0, "A":0,"a":0 ,"G":0, "g":0}
    numeric = []
    for base in sequence:
        numeric.append(dict_of_bases[base])
    dft = fft(np.array(numeric))
    mag = abs(dft)
    mag_avg = np.average(mag)
    return mag_avg


def magtropy(sequence):
    return magnitude_avg(sequence, representation = "Just A")/entropy(sequence)


# Saving Entropy values to dictionary
#removed temp_entropy_dict
file_path_1 = getcwd()
magtropy_dict = {}
for test in my_dict.keys():
    temp_dict = {}
    for family in  my_dict[test].keys():
        magtropy_values = []
        for file in my_dict[test][family]:
            start_seq = list(SeqIO.parse((f"{file_path_1}/data/{test}/{family}/{file}"), "fasta"))
            count = len(start_seq[0].seq)
            final_seq = "".join([char for char in start_seq[0].seq])
            magtropy_values.append((magtropy(final_seq)))
            temp_dict[family] = magtropy_values
    magtropy_dict[test] = temp_dict

magtropy_dict["Test1"]



colors = ["red", "royalblue", "gold", "darkorchid", "paleturquoise", "crimson", "m", "darkgreen", "olive", "aqua", "coral", "gray", "firebrick", "violet", "chartreuse"]
def plot_magtropy(test):
    loc = 0
    for family in magtropy_dict[test]:
        #sns.set(rc={"figure.figsize": (10, 7)})
        plt.figure(figsize = (8,6))
        print(sns.distplot(magtropy_dict[test][family], color=colors[loc], label=family, bins = 20))
        print(plt.title(family))
        loc = loc + 1
plot_magtropy("Test8")





    
