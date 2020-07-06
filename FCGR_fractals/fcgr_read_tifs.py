from Bio import SeqIO
import collections
from os import getcwd, listdir, makedirs, path, chdir
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import math
import json
from PIL import Image
from skimage.io import imread

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


# Reading files
k = 3
#for test in my_dict:
    #for folder in my_dict[test]:
for file in my_dict['Test3a']['Gammacoronavirus']:
    saved_path = getcwd() + f"/FCGR_fractals/plots/{k}-mers/Test3a/Gammacoronavirus"
    imread(saved_path + f'/{file[0:len(file)-6]}_{k}-mer_plot.tif')

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
                plt.title(str(k) + "-mer CGR:" + ' ' + file[0:len(file)-6])
                plt.imshow(chaos_kmer, interpolation='nearest', cmap=cm.gray_r)
                plt.show()
                saved_path = getcwd() + f"/FCGR_fractals/plots/{k}-mers/{test}/{folder}"
                plt.savefig(saved_path + f'/{file[0:len(file)-6]}_{k}-mer_plot.tif')
                #if not path.exists(saved_path):
                    #plt.savefig(saved_path + f'/{file[0:len(file)-6]}_{k}-mer_plot.tif')



saved_path = getcwd() + f"/FCGR_fractals/plots/7-mers/Test1/Anelloviridae"
imread(saved_path + '/Anelloviridae_199_7-mer_plot.tif')
