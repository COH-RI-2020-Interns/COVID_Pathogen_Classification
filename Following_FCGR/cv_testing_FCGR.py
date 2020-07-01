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
import collections
import cv2
import math
from matplotlib import cm
from PIL import Image
import pylab
from scipy.fftpack import dct, idct
from scipy.linalg import svd


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

def count_kmers(sequence, k):
    d = collections.defaultdict(int)
    for i in range(len(sequence)-(k-1)):
        d[sequence[i:i+k]] +=1
    for key in d.keys():
        if "N" in key:
            del d[key]
    return d


# getting the count of a specific kmer,
# dividing by (length of sequence - length of kmer + 1)
def probabilities(kmer_count, k, sequence):
        probabilities = collections.defaultdict(float)
        N = len(sequence)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - k + 1)
        return probabilities

def chaos_game_representation(probabilities, k):
        array_size = int(math.sqrt(4**k))
        chaos = []
        for i in range(array_size):
            chaos.append([0]*array_size)
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        for key, value in probabilities.items():
            for char in key:
                if char == "T":
                    posx += maxx / 2
                    posx = int(posx)
                elif char == "C":
                    posy += maxy / 2
                    posy = int(posy)
                elif char == "G":
                    posx += maxx / 2
                    posx = int(posx)
                    posy += maxy / 2
                    posy = int(posy)
                maxx = maxx / 2
                maxy /= 2
            chaos[posy-1][posx-1] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1

        return chaos




one = np.array([1,2,3,4,5])
dct = dct(one, type=3, norm='ortho')
dct


#___________________________________________________________________
# Printing PCC data by Test

new_dict_4 = {}

for folder in my_dict['Test3b']:
    list_sequences = []
    for file in my_dict["Test3b"][folder]:
        file_path = getcwd() + f"/data/Test3b/{folder}/{file}"
        seq  = make_sequence(file_path)
        f7 = count_kmers(seq, 7)
        f7_prob = probabilities(f7, 7, seq)
        chaos_k7 = chaos_game_representation(f7_prob, 7)
        x = np.array(chaos_k7)
        dct = dct(x[0], type=3, norm='ortho')
        U, s, VT = svd(dct)
        list_sequences.append(U, s, VT)
    new_dict_4[folder] = list_sequences

# run together for full plot
plt.title('Chaos game representation for 7-mers')
plt.imshow(chaos_k7, interpolation='nearest', cmap=cm.gray_r)
plt.show()


#___________________________________________________________________

#train_test_split
#Using delta and alphacoronavirus data from Test3b
delta = new_dict_4["Deltacoronavirus"]
alpha = new_dict_4["Alphacoronavirus"]
one = []
two  = []
for i in delta:
    one.append(len(i))
for i in alpha:
    two.append(len(i))
print(sorted(one))
print(sorted(two))

#making an array of all the names for the first column (Delta/Alpha)
name1 = []
for i in range(len(delta)):
    name1.append("1")
for i in range(len(alpha)-1):
    name1.append("2")


#dataframe has all the values of the magnitudes split up into base pairs
#inserting family name
df = pd.DataFrame(data= delta + alpha)
df = df.drop(68)
for i in df.columns:
    if(i>25401):
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
print(confusion_matrix(y_test,y_pred_k_neighbors))
print(classification_report(y_test,y_pred_k_neighbors))


#Linear SVM classifier:
linear_svm_classifier = SVC(kernel='linear')
linear_svm_classifier.fit(X_train, y_train)
y_pred_linear_svm = linear_svm_classifier.predict(X_test)
print(confusion_matrix(y_test,y_pred_linear_svm))
print(classification_report(y_test,y_pred_linear_svm))



#Linear Discriminant classifier:
linear_discriminant_classifier = LinearDiscriminantAnalysis()
linear_discriminant_classifier.fit(X_train, y_train)
y_pred_linear_discriminant = linear_discriminant_classifier.predict(X_test)
print(confusion_matrix(y_test,y_pred_linear_discriminant))
print(classification_report(y_test,y_pred_linear_discriminant))


#Polynomial SVM Classifier (Types of SVM = linear, poly, rbf, etc)
polynomial_svm_classifier = SVC(kernel = "poly")
polynomial_svm_classifier.fit(X_train, y_train)
y_pred_polynomial_svm = polynomial_svm_classifier.predict(X_test)
print(confusion_matrix(y_test,y_pred_polynomial_svm))
print(classification_report(y_test,y_pred_polynomial_svm))



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



#Making K-mers of sequences and trying jaccard_similarity

seq1 = 'ATGGACCAGATATAGGGAGAGCCAGGTAGGACA'
seq2 = 'ATGGACCAGATATTGGGAGAGCCGGGTAGGACA'


def k_mers(sequence, length):
    k_mers_seq = []
    for i in range(0,len(sequence)-length+1):
        k_mers_seq.append(sequence[i:i+3])
    return k_mers_seq

def jaccard_similarity(a, b):
    a = set(a)
    b = set(b)

    intersection = len(a.intersection(b))
    union = len(a.union(b))

    return intersection / union

K = 10
kmers1 = k_mers(seq1, K)
kmers2 = k_mers(seq2, K)
print(jaccard_similarity(kmers1, kmers2))
