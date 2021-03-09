import numpy as np
import pandas as pd
from os import getcwd, listdir, system

fd_path = getcwd() + f"/Work_Remaining/fd_magtropy_csvs"

csv_files = listdir(fd_path)
csv_files = sorted(csv_files)
csv_files

fd_df = pd.read_csv(f"{fd_path}/{csv_files[7]}")
fd_df = fd_df.drop(columns = ['Unnamed: 0', 'magtropy_PP'])
fd_df
fd_df['Sublevel Name'].value_counts()


magtropy_path = getcwd() + f"/Work_Remaining/all_feature_data_csvs"

csv_files2 = listdir(magtropy_path)
csv_files2 = sorted(csv_files2)
csv_files2

magtropy_df = pd.read_csv(f"{magtropy_path}/{csv_files2[7]}")
magtropy_df = magtropy_df.drop(columns = ['Unnamed: 0', 'Sublevel_Name'])
magtropy_df

full_df = pd.concat([fd_df, magtropy_df], axis=1)
full_df

# full_df1 = full_df[full_df['Sublevel Name'] == 'Cornidovirineae']
# full_df2 = full_df[full_df['Sublevel Name'] == 'Arnidovirineae']
# full_df2
# full_df = pd.concat([full_df1, full_df2])
# full_df

full_df['Sublevel Name'].value_counts()

taxonomic_level = input("taxonomic level: ")

saved_path = getcwd() + f"/Work_Remaining/all_feature_data_csvs"

full_df.to_csv(saved_path + f"/EIIP_{taxonomic_level[2:]}_100.csv")
