import numpy as np
from os import getcwd, listdir
from itertools import combinations

main_folder_path = getcwd() + "/data"

main_folders = sorted(listdir(main_folder_path))[1:3]

file_tuple_list = []

for folder in main_folders:
    for sub_folder in listdir(f"{main_folder_path}/{folder}"):
        for file in listdir(f"{main_folder_path}/{folder}/{sub_folder}"):
            file_tuple_list.append((folder, sub_folder, file))

file_combos = list(combinations(file_tuple_list, 2))

file_combos[0]
