import numpy as np
import json
from os import getcwd, listdir

file_path = getcwd() + "/data/JSON_Files"

f = open(f"{file_path}/{listdir(file_path)[0]}", )

my_dict = json.load(f)

my_dict 
