import numpy as np
import pandas as pd
from Bio import SeqIO
import rapidjson
from os import getcwd, listdir
from itertools import product
from scipy.fft import fft, ifft
from scipy import stats
import pywt
import matplotlib.pyplot as plt
import seaborn as sns
import json
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, classification_report

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
    rapidjson.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}", )

my_dict = rapidjson.load(f)


# Getting all possible combinations of 2 for the fasta files
file_tuple_list = []

for folder in folders:
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        for file in listdir(f"{folder_path}/{folder}/{sub_folder}"):
            file_tuple_list.append((folder, sub_folder,file))


file_combos = list(product(file_tuple_list, repeat=2))
new_dict = {}
#each of the keys shows tuples, first element is the virus, second is file
#new_dict contains all the fasta files in each family in each test
for i in my_dict.keys():
    file_list = []
    for j in my_dict[i].keys():
        for file in my_dict[i][j]:
            file_list.append((j,file))
        new_dict[i] = file_list

#new_dict2 contains all the products of the fasta files in the test folder
new_dict_2 = {}
for key in new_dict.keys():
    seq_perm = list(product(new_dict[key], repeat=2))
    new_dict_2[key] =  seq_perm


len(new_dict_2["Test6"])
48*48

len(new_dict_2["Test3a"])
82 * 82

len(new_dict_2["Test3b"])
74 * 74



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

# symmetric padding
def normalization(numerical1, numerical2):
    if len(numerical1)>len(numerical2):
        dna_seq = np.array(numerical2)
        numer = len(numerical1)- len(numerical2)
        if(numer%2 != 0):
            numer = numer - 0.5
            pad_width = numer/2
            numerical2 = pywt.pad(dna_seq, pad_width, "symmetric")
            numerical1 = numerical1[0:len(numerical1)-1]
        else:
            pad_width = numer/2
            numerical2 = pywt.pad(dna_seq,pad_width, "symmetric")
    elif len(numerical1)<len(numerical2):
        dna_seq = np.array(numerical1)
        numer = len(numerical2)- len(numerical1)
        if(numer%2 != 0):
            numer = numer - 0.5
            pad_width = numer/2
            numerical1 = pywt.pad(dna_seq, pad_width, "symmetric")
            numerical2 = numerical2[0:len(numerical2)-1]
        else:
            pad_width = numer/2
            numerical1 = pywt.pad(dna_seq,pad_width, "symmetric")
    else:
        numerical1 = numerical1
        numerical2 = numerical2

    return numerical1,numerical2
# CHECK WHY ANTISYMMETRIC PADDING DOESN"T WORK


#uses new_dict_2 and finds all the magnitudes of each fasta file in the tuples
#result  = np.array((mag_file1,mag_file2)), np.array((mag_file1, mag_file3))
def magnitude_array(test, dict):
    mag_list = []
    for tuple in dict[test]:
        file_path = getcwd() + f"/data/{test}/{tuple[0][0]}/{tuple[0][1]}"
        file_path2 = getcwd() + f"/data/{test}/{tuple[1][0]}/{tuple[1][1]}"
        seq1  = make_sequence(file_path)
        seq2 = make_sequence(file_path2)
        pp1 = numerical_pp(seq1)
        pp2 = numerical_pp(seq2)
        pp1_norm = normalization(pp1,pp2)[0]
        pp2_norm = normalization(pp1,pp2)[1]
        fft_1 = fft(pp1_norm)
        fft_2 = fft(pp2_norm)
        mag_1 = abs(fft_1)
        mag_2 = abs(fft_2)
        mag_list.append(np.array((mag_1, mag_2)))
    return mag_list

test3a = magnitude_array("Test3a", new_dict_2)
test6 = magnitude_array("Test6", new_dict_2)

#uses the two values to make one list of pearsons correlations
def pearsons(magnitude_array):
    pearson_corr = []
    for perm in magnitude_array:
        pearson_corr.append((stats.pearsonr(perm[0], perm[1]))[0])
    return pearson_corr

#Reshaping the array to be of size 82 by 82
test3a_pearsons = pearsons(test3a)
test3a_pearsons = np.array(test3a_pearsons).reshape(82,82)
test6_pearsons = pearsons(test6)
test6_pearsons = np.array(test3a_pearsons).reshape(48,48)



len(test3a_pearsons)
test3a_pearsons.shape

# Getting Clusters
data_path = getcwd() + "/data/JSON_Files"

with open(f"{data_path}/{listdir(data_path)[0]}", "r") as my_file:
    my_fasta = rapidjson.load(my_file)

my_fasta.keys()

def get_target(test_name):
    cluster_info = []
    cluster_name = []
    for key in my_fasta[test_name].keys():
        for file_name in my_fasta[test_name][key]:
             cluster_name.append(key)
             cluster_info.append(len(my_fasta[test_name][key]))
    return cluster_name

# Hypertuning
model_dict = {'log': LogisticRegression(),
             'rf': RandomForestClassifier(),
             'ada': AdaBoostClassifier(),
             'knn': KNeighborsClassifier(),
             'svm': SVC()
                }

data_path = getcwd() + "/data/JSON_Files"

#opening the json file that contains all the different parameters of each classification model
with open(f"{data_path}/{(listdir(data_path))[1]}", "r") as f:
    parameter_config = json.load(f)
parameter_config


def ML_Pipeline(features, target, estimator, cv, test_size, print_results=None):

    # Split Data into Training and Testing
    X_train, X_test, Y_train, Y_test = train_test_split(features, target, test_size=test_size, stratify=target)

    # Creating a Hyperparameter Tuning Strategy
    base_model = model_dict[estimator]
    model_params = parameter_config[estimator]
    ml_model = RandomizedSearchCV(base_model, model_params, cv=cv)

    # Train the Model
    ml_model.fit(X_train, Y_train)

    # Getting the Best Parameters and Results from Cross Validation
    print(f"The best algorithm is: {ml_model.best_estimator_} \n")
    print(f"The mean cross-validated score is: {ml_model.best_score_} \n")
    print(f"The best parameters for this model is: {ml_model.best_params_} \n")

    if print_results == 'yes':
        print(f"The cross validation results are: {ml_model.cv_results_}")

    # Getting Predictions on Holdout Set
    Y_pred = ml_model.predict(X_test)

    # Evaluating the Holdout
    print(f"The accuracy score is: {round(accuracy_score(Y_test, Y_pred)*100, 2)}% \n")
    print(f"The Matthew's Correlation Coefficient is: {matthews_corrcoef(Y_test, Y_pred)} \n")
    print(f"The confusion matrix is: {confusion_matrix(Y_test, Y_pred)} \n")
    print(classification_report(Y_test, Y_pred))

    return ml_model


#def ML_Pipeline(features, target, estimator, cv, test_size, print_results=None):
ML_Pipeline(test3a_pearsons, get_target("Test3a"), "knn" , 10, 0.2, print_results=True)




#why did antisymmetric padding not work
#why are we getting different result every time we run the ml pipeline
