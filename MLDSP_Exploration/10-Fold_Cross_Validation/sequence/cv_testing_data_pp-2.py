import numpy as np
import pandas as pd
import json, pywt, math
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn import datasets, svm, metrics
from os import getcwd, listdir, system
from itertools import permutations
from scipy.fft import fft, ifft
from scipy import stats
from Bio import SeqIO
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier


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

for folder in my_dict['Test1']:
    list_sequences = []
    for file in my_dict["Test1"][folder]:
        file_path = getcwd() + f"/data/Test1/{folder}/{file}"
        seq  = make_sequence(file_path)
        pp = numerical_pp(seq)
        list_sequences.append(pp)
    new_dict_4[folder] = list_sequences




#___________________________________________________________________

#train_test_split
#Using delta and alphacoronavirus data from Test3b
anelloviridae = new_dict_4["Anelloviridae"]
genomoviridae = new_dict_4["Genomoviridae"]
microviridae = new_dict_4["Microviridae"]
ortervirales = new_dict_4["Ortervirales"]
parvoviridae = new_dict_4["Parvoviridae"]

list_one= [anelloviridae, genomoviridae, microviridae, ortervirales, parvoviridae]
anelloviridae
two = []
for i in list_one:
    one = []
    for j in i:
        one.append(len(j))
    two.append(min(one))
two


#making an array of all the names for the first column (Delta/Alpha)
name1 = []
for i in range(len(anelloviridae)):
    name1.append("1")
for i in range(len(genomoviridae)):
    name1.append("2")
for i in range(len(microviridae)):
    name1.append("3")
for i in range(len(ortervirales)):
    name1.append("4")
for i in range(len(parvoviridae)):
    name1.append("5")

len(name1)

#dataframe has all the values of the magnitudes split up into base pairs
#inserting family name
df = pd.DataFrame(data= anelloviridae + genomoviridae + microviridae + ortervirales+ parvoviridae)
for i in df.columns:
    if(i>2001):
        df = df.drop(columns = [i])
df.insert(0, "Family", name1)
df

#Setting X to magnitudes, y to be family name
y = df["Family"]
X = df.drop(columns = ["Family"])
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, shuffle = True, random_state=0) # making Test size 0.2 instead of 0.1

#Standardizing the scale of the X train and X test to fall between -1 and 1
sc_X = StandardScaler()
X_train = sc_X.fit_transform(X_train)
X_test  = sc_X.transform(X_test)



#K_neighbors classification:
k_value = int(math.sqrt(len(y_test)) )#using a k value of 3, odd number and closest to
k_neighbors_classifier = KNeighborsClassifier(n_neighbors = k_value, p = 2, metric = "euclidean")
k_neighbors_classifier.fit(X_train, y_train)   #fitting the classifier on the training data, testing the ouput with y-pred
y_pred_k_neighbors = k_neighbors_classifier.predict(X_test)
y_pred_k_neighbors
print("CM:" , confusion_matrix(y_test,y_pred_k_neighbors))
print("CR:" , classification_report(y_test,y_pred_k_neighbors))
print("Accuracy score:" ,  accuracy_score(y_test,y_pred_k_neighbors))


#Linear SVM classifier:
linear_svm_classifier = SVC(kernel='linear')
linear_svm_classifier.fit(X_train, y_train)
y_pred_linear_svm = linear_svm_classifier.predict(X_test)
print("CM:" , confusion_matrix(y_test,y_pred_linear_svm))
print("CR:" , classification_report(y_test,y_pred_linear_svm))
print("Accuracy score:" ,  accuracy_score(y_test,y_pred_linear_svm))



#Linear Discriminant classifier:
linear_discriminant_classifier = LinearDiscriminantAnalysis()
linear_discriminant_classifier.fit(X_train, y_train)
y_pred_linear_discriminant = linear_discriminant_classifier.predict(X_test)
print("CM:" , confusion_matrix(y_test,y_pred_linear_discriminant))
print("CR:" , classification_report(y_test,y_pred_linear_discriminant))
print("Accuracy score:" ,  accuracy_score(y_test,y_pred_linear_discriminant))


#Polynomial SVM Classifier (Types of SVM = linear, poly, rbf, etc)
polynomial_svm_classifier = SVC(kernel = "poly")
polynomial_svm_classifier.fit(X_train, y_train)
y_pred_polynomial_svm = polynomial_svm_classifier.predict(X_test)
y_pred_polynomial_svm
print("CM:" , confusion_matrix(y_test,y_pred_polynomial_svm))
print("CR:" , classification_report(y_test,y_pred_polynomial_svm))
print("Accuracy score:" ,  accuracy_score(y_test,y_pred_polynomial_svm))

#Random Forest classifier
random_forest_classifier = RandomForestClassifier(n_estimators = 100)
random_forest_classifier.fit(X_train, y_train)
y_pred_random_forest = random_forest_classifier.predict(X_test)
print("CM:" , confusion_matrix(y_test,y_pred_random_forest))
print("CR:" , classification_report(y_test,y_pred_random_forest))
print("Accuracy score:" ,  accuracy_score(y_test,y_pred_random_forest))

feature_imp= pd.DataFrame({'Value':model.feature_importance(),'Feature':X.columns})
plt.figure(figsize=(40, 20))
sns.set(font_scale = 5)
sns.barplot(x="Value", y="Feature", data=feature_imp.sort_values(by="Value", ascending=False)[0:20])
plt.title('LightGBM Features (avg over folds)')
plt.tight_layout()
plt.savefig('lgbm_importances-01.png')
plt.show()


# Basic shaping of the linear and quadratic svm graphs
X = np.c_[(.4, -.7),
          (-1.5, -1),
          (-1.4, -.9),
          (-1.3, -1.2),
          (-1.1, -.2),
          (-1.2, -.4),
          (-.5, 1.2),
          (-1.5, 2.1),
          (1, 1),
          # --
          (1.3, .8),
          (1.2, .5),
          (.2, -2),
          (.5, -2.4),
          (.2, -2.3),
          (0, -2.7),
          (1.3, 2.1)].T
Y = [0] * 8 + [1] * 8



# figure number
fignum = 1

# fit the model
for kernel in ('linear', 'poly'):
    clf = svm.SVC(kernel=kernel, gamma=2)
    clf.fit(X, Y)

    # plot the line, the points, and the nearest vectors to the plane
    plt.figure(fignum, figsize=(4, 3))
    plt.clf()

    plt.scatter(clf.support_vectors_[:, 0], clf.support_vectors_[:, 1], s=80,
                facecolors='none', zorder=10, edgecolors='k')
    plt.scatter(X[:, 0], X[:, 1], c=Y, zorder=10, cmap=plt.cm.Paired,
                edgecolors='k')

    plt.axis('tight')
    x_min = -3
    x_max = 3
    y_min = -3
    y_max = 3

    XX, YY = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
    Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(XX.shape)
    plt.figure(fignum, figsize=(4, 3))
    plt.pcolormesh(XX, YY, Z > 0, cmap=plt.cm.Paired)
    plt.contour(XX, YY, Z, colors=['k', 'k', 'k'], linestyles=['--', '-', '--'],
                levels=[-.5, 0, .5])

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    plt.xticks(())
    plt.yticks(())
    fignum = fignum + 1
plt.show()
