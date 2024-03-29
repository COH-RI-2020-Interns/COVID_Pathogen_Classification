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
from sklearn.tree import DecisionTreeClassifier


#Going to Test folders
folder_path = getcwd() + "/data2"

folders = sorted(listdir(folder_path))[0:11]
folders

folder_dict = {}

# Going to Specific Virus Folders inside the Test folders
for folder in folders:
    subfolder_dict = {}
    for sub_folder in listdir(f"{folder_path}/{folder}"):
        subfolder_dict[sub_folder] = listdir(f"{folder_path}/{folder}/{sub_folder}")
    folder_dict[folder] = subfolder_dict


#Adding data to a JSON file
output_path = getcwd() + "/data2/JSON_Files"

with open(f"{output_path}/fasta_files.json", "w") as my_file:
    json.dump(folder_dict, my_file)

f = open(f"{output_path}/{listdir(output_path)[0]}")

my_dict = json.load(f)
for test in my_dict:
    test = sorted(test)

#Dictionary of numerical representations
rep_dict = {"Int1":{"T":0,"t":0,"C":1,"c":1, "A":2,"a":2 ,"G":3, "g":3},
"Int2": {"T":1,"t":1,"C":2,"c":2, "A":3,"a":3 ,"G":4, "g":4},
"Real": {"T":-1.5,"t":-1.5,"C":0.5,"c":0.5, "A":1.5,"a":1.5 ,"G":-1.5, "g":-1.5},
"EIIP": {"T":0.1335,"t":0.1335,"C":0.1340,"c":0.1340, "A":0.1260,"a":0.1260 ,"G":0.0806, "g":0.0806},
"PP": {"T":1,"t":1,"C":1,"c":1, "A":-1,"a":-1 ,"G":-1, "g":-1},
"Paired Numeric": {"T":1,"t":1,"C":-1,"c":-1, "A":1,"a":1 ,"G":-1, "g":-1},
"Just A": {"T":0,"t":0,"C":0,"c":0, "A":1,"a":1 ,"G":0, "g":0},
"Just C": {"T":0,"t":0,"C":1,"c":1, "A":0,"a":0 ,"G":0, "g":0},
"Just G": {"T":0,"t":0,"C":0,"c":0, "A":0,"a":0 ,"G":1, "g":1},
"Just T": {"T":1,"t":1,"C":0,"c":0, "A":0,"a":0 ,"G":0, "g":0}}


my_dict["1_Realm"]
# Finding the Average Magnitude of the Sequence
def magnitude_avg(sequence):
    mag_avg_list = []
    base_list = ["D", "K", "M", "N", "R", "S", "W", "Y"]
    for i in base_list:
        sequence = sequence.replace(i, "")
    for rep in rep_dict:
        dict_of_bases = rep_dict[rep]
        numeric = []
        for base in sequence:
            numeric.append(dict_of_bases[base])
        dft = fft(np.array(numeric))
        mag = abs(dft)
        mag_avg = np.average(mag)
        mag_avg_list.append(mag_avg)
    return mag_avg_list

# Calculating Entropy
def entropy(sequence):
    counts = Counter(sequence)
    props = {key: counts[key] / sum(counts.values()) for key in counts}
    products = {key: props[key]*np.log(props[key]) for key in props}
    return -1 * sum(products.values())

#Calculating Magtropy
def magtropy(sequence):
    list_magtropy = [avg/entropy(sequence) for avg in magnitude_avg(sequence)]
    return list_magtropy


# Saving Magtropy values to dictionary for specific sublevel
def magtropy_dict(sublevel):
    file_path_1 = getcwd()
    magtropy_dict = {}
    magtropy_values = []
    for folder in my_dict[sublevel].keys():
        for file in my_dict[sublevel][folder]:
            start_seq = list(SeqIO.parse((f"{file_path_1}/data2/{sublevel}/{folder}/{file}"), "fasta"))
            final_seq = "".join([char for char in start_seq[0].seq])
            magtropy_values.append((folder, magtropy(final_seq)[0]))
            #If you would like to use more representations, you can add it in with magtropy(final_seq)[index]
            #If you would like to change the representation, change the index number, refer to rep. dict
    magtropy_dict[sublevel] = magtropy_values
    sublevel = pd.DataFrame.from_dict (magtropy_dict[sublevel])
    sublevel.columns = ["Sublevel Name", "Int1", "Int2", "Real", "EIIP", "PP", "Paired Numeric", "JustA", "JustC", "JustG", "JustT"]
    return sublevel




#Preparing training data for supervised machine learning
sublevel_df = magtropy_dict("10_Subgenus")

#sublevel_df


X = sublevel_df.drop(columns = ["Sublevel Name"])    #these are the training features
y = pd.DataFrame(sublevel_df["Sublevel Name"])       #this are the target labels

sublevel_df
# Hypertuning
model_dict = {'log': LogisticRegression(),
             'rf': RandomForestClassifier(),
             'ada': AdaBoostClassifier(),
             'knn': KNeighborsClassifier(),
             'svm': SVC(),
             'decision_tree': DecisionTreeClassifier()
                }


data_path = getcwd() + "/data2/JSON_Files"


#opening the json file that contains all the different parameters of each classification model
with open(f"{data_path}/{(listdir(data_path))[1]}", "r") as f:
    parameter_config = json.load(f)

def ML_Pipeline(features, target, estimator, cv, test_size, print_results=None):

    # Split Data into Training and Testing
    X_train, X_test, Y_train, Y_test = train_test_split(features, target, test_size=test_size, stratify=target)

    # Creating a Hyperparameter Tuning Strategy
    if estimator == "svm":
        model_dict[estimator].probability = True
    base_model = model_dict[estimator]
    model_params = parameter_config[estimator]
    ml_model = RandomizedSearchCV(base_model, model_params, n_iter= 15, cv=cv)

    # Train the Model
    ml_model.fit(X_train, np.ravel(Y_train))

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


my_model = ML_Pipeline(X, y, "svm", 10, 0.2, print_results = 'yes')



#Testing data of COVID-19 Files
covid_df = magtropy_dict("0_COVID")
X_test = covid_df.drop(columns = ["Sublevel Name"])    #these are the testing features
my_model.predict(X_test)



X_test = covid_df.drop(columns = ["Sublevel Name"])    #these are the testing features
my_model.predict(X_test)
# getting 1file COVID sequences


covid_df2 = pd.read_csv(getcwd() + "/covid.csv")
X_test2 = covid_df2.drop(columns = ["Sublevel Name"])
my_model.predict(X_test2)
#
#
# classes = my_model.classes_
# classes
#
# probas = my_model.predict_proba(X_test)
# [(i, np.where(probas == i)) for i in probas if (i[3]>i[2] and i[0]<i[3])]
# [(i, np.where(probas == i)) for i in probas if i[3]>i[2]]

#my_model.predict_proba(X1)
