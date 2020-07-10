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

<<<<<<< HEAD
folders = sorted(listdir(folder_path))[1:9]
=======
folders = sorted(listdir(folder_path))[2:9]
>>>>>>> b4ea0f3bcf166a147c97004301490374b69ebe33
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

my_dict
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
    return magnitude_avg(sequence)/entropy(sequence)


# Saving Entropy values to dictionary
#removed temp_entropy_dict
file_path_1 = getcwd()
entropy_dict = {}
for test in my_dict.keys():
    entropy_values = []
    for family in  my_dict[test].keys():
        for file in my_dict[test][family]:
            start_seq = list(SeqIO.parse((f"{file_path_1}/data/{test}/{family}/{file}"), "fasta"))
            count = len(start_seq[0].seq)
            final_seq = "".join([char for char in start_seq[0].seq])
            entropy_values.append((family,magtropy(final_seq)))
    entropy_dict[test] = entropy_values
entropy_dict["Test7"]
test1 = pd.DataFrame.from_dict(entropy_dict["Test1"])
test2 = pd.DataFrame.from_dict(entropy_dict["Test2"])
test3a = pd.DataFrame.from_dict(entropy_dict["Test3a"])
test3b = pd.DataFrame.from_dict(entropy_dict["Test3b"])
test4 = pd.DataFrame.from_dict(entropy_dict["Test4"])
test5  = pd.DataFrame.from_dict(entropy_dict["Test5"])
test6 = pd.DataFrame.from_dict(entropy_dict["Test6"])
test7 = pd.DataFrame.from_dict(entropy_dict["Test7"])

test1.columns = ["Family", "Magtropy"]
test2.columns = ["Family", "Magtropy"]
test3a.columns = ["Family", "Magtropy"]
test3b.columns = ["Family", "Magtropy"]
test4.columns = ["Family", "Magtropy"]
test5.columns = ["Family", "Magtropy"]
test6.columns = ["Family", "Magtropy"]
test7.columns = ["Family", "Magtropy"]

# Hypertuning
model_dict = {'log': LogisticRegression(),
             'rf': RandomForestClassifier(),
             'ada': AdaBoostClassifier(),
             'knn': KNeighborsClassifier(),
             'svm': SVC()
                }
def removeCovid(test):
    test = test.drop([test[test["Family"] == "COVID19"].index[0]], axis = 0)
    return test


test5 = removeCovid(test5)
test5[test5["Family"] == "COVID19"]

#df=pd.DataFrame(test4["Family"])

X = pd.DataFrame(test7["Magtropy"])
X
y = pd.DataFrame(test7["Family"])
y
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
    ml_model = RandomizedSearchCV(base_model, model_params, n_iter= 15, cv=cv)

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


my_model = ML_Pipeline(X, y, "knn", 10, 0.2, print_results = None)



start_seq = list(SeqIO.parse((getcwd() + f"/data/Test5/COVID19/MN908947.fasta"), "fasta"))
count = len(start_seq[0].seq)
final_seq = "".join([char for char in start_seq[0].seq])
COVID = magtropy(final_seq)

COVID

COVID = {"Magtropy": [COVID]}
COVID
df2 = pd.DataFrame(COVID, columns = ["Magtropy"])
df2

my_model.predict(df2)
