from Bio import SeqIO
import collections
from os import getcwd, listdir, makedirs, path, chdir
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import math
import json

# Accessing the Fasta Files
folder_path = getcwd() + "/data"

folders = sorted(listdir(folder_path))[1:8]
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


# Extracting the Sequence from the Fasta File
def make_sequence(path_of_file):
    start_seq = list(SeqIO.parse(f"{path_of_file}", "fasta"))
    count = len(start_seq[0].seq)
    final_seq = "".join([char for char in start_seq[0].seq])
    return final_seq

# Counting total k-mers possible
def count_kmers(final_seq, k):
    d = collections.defaultdict(int)
    for i in range(len(final_seq)-(k-1)):
        d[final_seq[i:i+k]] +=1
    for key in d.keys():
        if "N" in key:
            del d[key]
    return d
#The N is used to remove a certain base pair A,T,C,G if needed
#finding the number of each unique k-mer in the sequence, for example
#it finds the number of time the kmer CAT appears in the sequence, and continues with each 3-mer


# getting the count of a specific kmer, dividing by (length of sequence - length of kmer +1)
def probabilities(kmer_count, k, final_seq):
        probabilities = collections.defaultdict(float)
        N = len(final_seq)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - k + 1)
        return probabilities

# Finding the Frequency of kmers
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

# Creating Folders to Save Plots
for test in my_dict:
    for folder in my_dict[test]:
        for k in range(3,8):
            saved_path = getcwd() + f"/FCGR_fractals/plots/{k}-mers/{test}/{folder}"
            if not path.exists(saved_path):
                makedirs(saved_path)

# Generating the Chaos Game Representation Plots
for k in range(3,8):
    for test in my_dict:
        for folder in my_dict[test]:
            for file in my_dict[test][folder]:
                file_path = getcwd() + f"/data/{test}/{folder}/{file}"
                seq  = make_sequence(file_path)
                freq = count_kmers(seq, k)
                freq_prob = probabilities(freq, k, seq)
                chaos_kmer = chaos_game_representation(freq_prob, k)
                #plt.title(str(k) + "-mer CGR:" + ' ' + file[0:len(file)-6])
                plt.imshow(chaos_kmer, interpolation='nearest', cmap=cm.gray_r)
                plt.show()
                saved_path = getcwd() + f"/FCGR_fractals/plots/{k}-mers/{test}/{folder}"
                plt.savefig(saved_path + f'/{file[0:len(file)-6]}_{k}-mer_plot.tif')
                # if not path.exists(saved_path):
                #     plt.savefig(saved_path + f'/{file[0:len(file)-6]}_{k}-mer_plot.tif')
