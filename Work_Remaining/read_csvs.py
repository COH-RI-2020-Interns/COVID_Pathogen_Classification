
import numpy as np
import pandas as pd
from os import getcwd, listdir, system

path = getcwd() + f"/Work_Remaining/fd_magtropy_csvs"

csv_files = listdir(path)
csv_files = sorted(csv_files)
csv_files

df = pd.read_csv(f"{path}/{csv_files[3]}")
df
df['Sublevel Name'].value_counts()
