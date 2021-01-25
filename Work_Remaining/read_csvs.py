
import numpy as np
import pandas as pd
from os import getcwd, listdir, system

path = getcwd() + f"/fd_magtropy_csvs"

csv_files = listdir(path)
csv_files = sorted(csv_files[:-1])

df = pd.read_csv(f"{path}/{csv_files[0]}")
df
df['Sublevel Name'].value_counts()
