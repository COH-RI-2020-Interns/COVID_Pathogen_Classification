import numpy as np
import pandas as pd
from collections import Counter
from os import getcwd, listdir



test = input("Please input the test you would like data from")
folder = input("Please input the folder you would like data from")

file_path = getcwd() + "/COVID_Pathogen_Classification/data"
