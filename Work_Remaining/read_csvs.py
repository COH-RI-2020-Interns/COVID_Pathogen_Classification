
import numpy as np
import pandas as pd
from os import getcwd, listdir, system

path = getcwd() + f"/Work_Remaining/embecovirus_EIIP_results"

csv_files = listdir(path)
csv_files = sorted(csv_files)
csv_files

df = pd.read_csv(f"{path}/{csv_files[13]}")

df
