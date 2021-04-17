
import numpy as np
import pandas as pd
from os import getcwd, listdir, system

taxonomic_level = input("Taxonomic Level: ")
path = getcwd() + f"/Work_Remaining/confidence_scores/{taxonomic_level}"

csv_files = listdir(path)
csv_files = sorted(csv_files)
csv_files

magtropy_df = pd.read_csv(f"{path}/{csv_files[2]}")
magtropy_df = magtropy_df.drop(columns='Unnamed: 0')
magtropy_df.to_csv(path + f"/magtropy_{taxonomic_level}.csv")

test_magtropy = magtropy_df[magtropy_df['Sublevel Name'] != magtropy_df['Label']]
len(test_magtropy.index)



fd_df = pd.read_csv(f"{path}/{csv_files[1]}")
fd_df = fd_df.drop(columns='Unnamed: 0')
fd_df.to_csv(path + f"/fd_{taxonomic_level}.csv")

test_fd = fd_df[fd_df['Sublevel Name'] != fd_df['Label']]
len(test_fd.index)

fd_magtropy_df = pd.read_csv(f"{path}/{csv_files[0]}")
fd_magtropy_df = fd_magtropy_df.drop(columns='Unnamed: 0')
fd_magtropy_df.to_csv(path + f"/fd+magtropy_{taxonomic_level}.csv")

test_fd_magtropy = fd_magtropy_df[fd_magtropy_df['Sublevel Name'] != fd_magtropy_df['Label']]
len(test_fd_magtropy.index)
